#pragma once

#include <cinttypes>
#include <thread>
#include <vector>
#include <chrono>
#include <atomic>
#include <exception>
#include <mutex>
#include <array>
#include <functional>
#include <string>
#include <future>

#include <iostream>

//#define REFRESH_ATP_DEBUG
//#define REFRESH_ATP_STATS

namespace refresh 
{
	namespace utils
	{
		// *********************************************************************
		inline void noop()
		{
#if defined(_WIN32)
			_mm_pause();
#elif defined(__aarch64__)
			std::this_thread::yield();
#else
			__builtin_ia32_pause();
#endif
		}

		// *********************************************************************
		class spin_mutex
		{
			std::atomic_flag a_lock{};

		public:
			// *********************************************************************
			void lock()
			{
				while (a_lock.test_and_set(std::memory_order_acquire))
					while (a_lock.test(std::memory_order_relaxed))
						;
			}

			// *********************************************************************
			void lock_with_noop()
			{
				while (a_lock.test_and_set(std::memory_order_acquire))
					while (a_lock.test(std::memory_order_relaxed))
						utils::noop();
			}

			// *********************************************************************
			void unlock()
			{
				a_lock.clear();
			}
		};
	}

	// *********************************************************************
	class active_thread_pool_state
	{
		friend class active_thread_pool;

		std::atomic_flag is_operating{};
		std::atomic<size_t> no_working{};
		std::exception_ptr exception_ptr{};

		std::mutex mtx;

		void inc()
		{
			no_working.fetch_add(1);
			is_operating.test_and_set();
		}

		void dec()
		{
			if (no_working.fetch_sub(1) == 1)
				is_operating.clear();

			is_operating.notify_all();
		}

		void set_exception(std::exception_ptr ptr)
		{
			std::lock_guard<std::mutex> lck(mtx);
			exception_ptr = ptr;
		}

	public:
		active_thread_pool_state() = default;
		
		void reserve(size_t n)	
		{
			// Do nothing here
		}
			
		std::exception_ptr get_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr;
		}

		bool was_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr != nullptr;
		}

		void wait()
		{
			is_operating.wait(true);
		}

		void heavy_wait()
		{
			while (is_operating.test())
				;
		}

		void busy_wait()
		{
			while (is_operating.test())
				refresh::utils::noop();
		}

		bool is_ready()
		{
			return !is_operating.test();
		}
	};

	// *********************************************************************
	class active_thread_pool
	{
	public:
		using pool_state_t = active_thread_pool_state;

	private:
		size_t init_pool_size;
		size_t max_pool_size;
		size_t busy_ticks;
		size_t no_workers{};
		size_t prev_worker_id{};
		std::chrono::duration<double> busy_time;

		utils::spin_mutex mtx;

		static const uint64_t thread_state_sleeping = 0;
		static const uint64_t thread_state_waiting = 1;
		static const uint64_t thread_state_code_loading = 2;
		static const uint64_t thread_state_working = 3;
		static const uint64_t thread_state_terminating = 4;

#ifdef REFRESH_ATP_STATS
		std::atomic<uint64_t> n_sleeping_use{};
		std::atomic<uint64_t> n_waiting_use{};
		std::atomic<uint64_t> n_waiting_to_sleeping{};
#endif

		// *********************************************************************
		struct worker_t
		{
			alignas(uint64_t) uint64_t state;
			std::thread internal_thread;
			size_t internal_id;
			std::function<void()> task{};
			pool_state_t* pool_ptr;
			size_t busy_ticks{};
			std::array<uint8_t, 64> fill{};				// to avoid false sharing

/*#ifdef REFRESH_ATP_STATS
			std::atomic<uint64_t> loc_sleeping_use{};
			std::atomic<uint64_t> loc_waiting_use{};
			std::atomic<uint64_t> loc_waiting_to_sleeping{};
#endif*/

		private:
			// *********************************************************************
			void worker_body()
			{
				std::atomic_ref<uint64_t> thread_state(state);

				while (true)
				{
					auto current_state = thread_state.load();

#ifdef REFRESH_ATP_DEBUG
					std::cerr << "my state: " + std::to_string(internal_id) + " - " + std::to_string(current_state) + "\n";
#endif

					switch (current_state)
					{
					case thread_state_sleeping:
						thread_state.wait(thread_state_sleeping);
						break;
					case thread_state_terminating:
						return;
					case thread_state_code_loading:
						while (thread_state.load() == thread_state_code_loading)
							utils::noop();
						break;
					case thread_state_waiting:
						for (size_t i = 0; i < busy_ticks; ++i)
						{
							if (thread_state.load() != thread_state_waiting)
								break;
							utils::noop();
						}
						
						{
							auto tmp = thread_state_waiting;
							thread_state.compare_exchange_strong(tmp, thread_state_sleeping);

#ifdef REFRESH_ATP_STATS
#else
#endif
						}
						break;
					case thread_state_working:
						try
						{
							task();

							if (pool_ptr)
								pool_ptr->dec();
						}
						catch (const std::exception &e)
						{
							if (pool_ptr)
							{
								pool_ptr->dec();

								std::cerr << "Exception: " + std::string(e.what()) + "\n";
								fflush(stderr);

								pool_ptr->set_exception(std::current_exception());
								
								return;
							}
						}

						{
//							auto tmp = thread_state_working;
//							thread_state.compare_exchange_strong(tmp, thread_state_waiting);
							thread_state.store(thread_state_waiting);
						}
						break;
					}
				}
			}

		public:
			// *********************************************************************
			worker_t() = delete;

			worker_t(size_t internal_id, size_t busy_ticks) :
				state(thread_state_waiting),
				internal_id(internal_id),
				pool_ptr(nullptr),
				busy_ticks(busy_ticks)
			{
				internal_thread = std::thread([&] {worker_body(); });
			}
		};

		// *********************************************************************
		std::vector<worker_t> workers;

		// *********************************************************************
		void empty_loop(size_t n)
		{
			volatile size_t aux = 0;

			for (size_t i = 0; i < n; ++i)
			{
				aux = aux + 1;
				refresh::utils::noop();
			}
		}

		// *********************************************************************
		void adjust_busy_time()
		{
			busy_ticks = 8;

			while (busy_ticks < (1 << 30))
			{
				auto t1 = std::chrono::high_resolution_clock::now();
				empty_loop(busy_ticks);
				auto t2 = std::chrono::high_resolution_clock::now();

				const std::chrono::duration<double> diff = t2 - t1;

				if (diff >= busy_time)
					break;

				if (busy_ticks & (busy_ticks - 1))
				{
					busy_ticks &= busy_ticks - 1;
					busy_ticks *= 2;
				}
				else
					busy_ticks += busy_ticks / 2;
			}
		}
		
		// *********************************************************************
		void init_workers()
		{
			workers.reserve(max_pool_size);

			for (size_t i = 0; i < init_pool_size; ++i)
				workers.emplace_back(i, busy_ticks);
		}

	public:
		// *********************************************************************
		active_thread_pool(size_t init_pool_size, size_t max_pool_size, std::chrono::duration<double> busy_time = std::chrono::milliseconds(2)) :
			init_pool_size(init_pool_size),
			max_pool_size(max_pool_size),
			busy_time(busy_time)
		{
			adjust_busy_time();

#ifdef REFRESH_ATP_DEBUG
			std::cout << "atp busy_ticks: " << busy_ticks << std::endl;
#endif

			init_workers();
		}

		// *********************************************************************
		active_thread_pool(size_t init_pool_size, size_t max_pool_size, size_t busy_ticks) :
			init_pool_size(init_pool_size),
			max_pool_size(max_pool_size),
			busy_ticks(busy_ticks)
		{
#ifdef REFRESH_ATP_DEBUG
			std::cout << "atp busy_ticks: " << busy_ticks << std::endl;
#endif

			init_workers();
		}

		// *********************************************************************
		~active_thread_pool()
		{
			cancel();
			for (auto& worker : workers)
				worker.internal_thread.join();
		}

		// *********************************************************************
		bool is_active()	const
		{
			return true;
		}

		// *********************************************************************
		void launch(std::function<void()>&& fun, active_thread_pool::pool_state_t* pool_state)
//		void launch(const std::function<void()>& fun, active_thread_pool::pool_state_t* pool_state)
		{
			mtx.lock();

			size_t n_workers = workers.size();
			size_t worker_id = (prev_worker_id + 1) % n_workers;
			size_t i = 0;

			// Find any waiting thread
			for (i = 0; i < n_workers; ++i)
			{
				uint64_t curr_state = thread_state_waiting;
				std::atomic_ref<uint64_t> ar(workers[worker_id].state);

				if (ar.compare_exchange_strong(curr_state, thread_state_code_loading))
				{
#ifdef REFRESH_ATP_STATS
					n_waiting_use++;
#endif
					break;
				}

				worker_id = (worker_id + 1 == n_workers) ? 0 : (worker_id + 1);
			}

			// If necessary look for sleeping thread
			if (i == n_workers)
			{
				for (i = 0; i < n_workers; ++i)
				{
					uint64_t curr_state = thread_state_sleeping;
					std::atomic_ref<uint64_t> ar(workers[worker_id].state);

					if (ar.compare_exchange_strong(curr_state, thread_state_code_loading))
					{
#ifdef REFRESH_ATP_STATS
						n_sleeping_use++;
#endif
						break;
					}

					worker_id = (worker_id + 1 == n_workers) ? 0 : (worker_id + 1);
				}
			}

#ifdef REFRESH_ATP_DEBUG
			std::cerr << "worker_id: " + std::to_string(worker_id) + "\n";
#endif

			// If necessary add new internal thread
			if (i == n_workers)
			{
				if (n_workers == max_pool_size)
				{
					mtx.unlock();
					// throw ...	
					exit(1);
				}

				worker_id = n_workers;

				workers.emplace_back(worker_id, busy_ticks);

				std::atomic_ref<uint64_t> ar(workers[worker_id].state);
				ar.store(thread_state_code_loading);
			}

			pool_state->inc();
			workers[worker_id].pool_ptr = pool_state;
			workers[worker_id].task = std::forward<std::function<void()>>(fun);
//			workers[worker_id].task = fun;

			std::atomic_ref<uint64_t> thread_state(workers[worker_id].state);

			thread_state.store(thread_state_working);
			thread_state.notify_one();

			prev_worker_id = worker_id;

			mtx.unlock();
		}

		// *********************************************************************
		void cancel()
		{
			mtx.lock();

			// Set state to terminating for all workers
			for (size_t i = 0; i < workers.size(); ++i)
			{
				std::atomic_ref<uint64_t> thread_state(workers[i].state);
				thread_state.store(thread_state_terminating);
				thread_state.notify_one();
				workers[i].internal_thread.join();

/*#ifdef REFRESH_ATP_STATS
				n_waiting_to_sleeping += workers[i].loc_waiting_to_sleeping;
#endif*/
			}

			workers.clear();

#ifdef REFRESH_ATP_STATS
			std::cerr << "No. waiting use        : " << n_waiting_use << std::endl;
			std::cerr << "No. sleeping use       : " << n_sleeping_use << std::endl;
			std::cerr << "No. waiting to sleeping: " << n_waiting_to_sleeping << std::endl;
#endif

			mtx.unlock();
		}
	};

	// *********************************************************************
	class async_pool_state
	{
		friend class async_pool;

		std::exception_ptr exception_ptr{};

		std::mutex mtx;

		std::vector<std::future<void>> workers;

		void add_future(std::future<void>&& fut)
		{
			std::lock_guard<std::mutex> lck(mtx);
			workers.emplace_back(std::move(fut));
		}

		void set_exception(std::exception_ptr ptr)
		{
			std::lock_guard<std::mutex> lck(mtx);
			exception_ptr = ptr;
		}

	public:
		async_pool_state() = default;

		void reserve(size_t n)
		{
			std::lock_guard<std::mutex> lck(mtx);
			workers.reserve(n);						
		}

		std::exception_ptr get_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr;
		}

		bool was_exception()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return exception_ptr != nullptr;
		}

		void wait()
		{
			std::lock_guard<std::mutex> lck(mtx);
			for (auto& f : workers)
			{
				try {
					f.wait();
				}
				catch (...)
				{
					set_exception(std::current_exception());
				}
			}

			workers.clear();
		}

		void heavy_wait()
		{
			wait();
		}

		void busy_wait()
		{
			wait();
		}

		bool is_ready()
		{
			std::lock_guard<std::mutex> lck(mtx);
			return workers.empty();
		}
	};

	// *********************************************************************
	class async_pool
	{
	public:
		using pool_state_t = async_pool_state;

	private:

	public:
		// *********************************************************************
		async_pool()
		{
		}

		// *********************************************************************
		~async_pool()
		{
		}

		// *********************************************************************
		bool is_active()	const
		{
			return true;
		}

		// *********************************************************************
		void launch(std::function<void()>&& fun, async_pool::pool_state_t* pool_state)
		{
			pool_state->add_future(std::async(std::launch::async, fun));
		}

		// *********************************************************************
		void cancel()
		{
			// Do nothing
		}
	};
}
