#pragma once
/*
This file is a part of Kmer-db software distributed under GNU GPL 3 licence.
The homepage of the Kmer-db project is http://sun.aei.polsl.pl/REFRESH/kmer-db

Authors: Sebastian Deorowicz, Adam Gudys, Maciej Dlugosz, Marek Kokot, Agnieszka Danek

*/

// Generic multithreading queues

#include <queue>
#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <vector>
#include <functional>
#include <cinttypes>

using namespace std;

#define USE_CPP20

#ifdef USE_CPP20
#include <atomic>

class Semaphore {
protected:
	atomic<int> counter;
	atomic_flag is_zero;

public:
	// *****************************************************************************************
	//
	Semaphore() : counter(0) {
		is_zero.clear();
	}

	Semaphore(const Semaphore& r) = delete;
	Semaphore& operator=(const Semaphore& r) = delete;

	// *****************************************************************************************
	//
	void inc() {
		is_zero.clear();
		counter.fetch_add(1);
	}

	// *****************************************************************************************
	//
	void inc(int num) {
		is_zero.clear();
		counter.fetch_add(num);
	}

	// *****************************************************************************************
	//
	void dec() {
		if (counter.fetch_sub(1) == 1)
		{
			is_zero.test_and_set();
			is_zero.notify_one();
		}
	}

	// *****************************************************************************************
	//
	void dec_notify_all() {
		if (counter.fetch_sub(1) == 1)
		{
			is_zero.test_and_set();
			is_zero.notify_all();
		}
	}

	// *****************************************************************************************
	//
	void waitForZero() {
		is_zero.wait(false);
	}

	void waitForZeroBusy()
	{
		while (!is_zero.test())
		{
#if defined(_WIN32)
			_mm_pause();
#elif defined(__aarch64__)
			std::this_thread::yield();
#else
			__builtin_ia32_pause();
#endif
		}
	}
};

#else
// *****************************************************************************************
//
class Semaphore {
protected:
	int counter;
	std::mutex mutex;
	std::condition_variable cv;

public:
	// *****************************************************************************************
	//
	Semaphore() : counter(0) {}

	// *****************************************************************************************
	//
	void inc() {
		std::unique_lock<std::mutex> lk(mutex);
		++counter;
	}

	// *****************************************************************************************
	//
	void inc(int num) {
		std::unique_lock<std::mutex> lk(mutex);
		counter += num;
	}

	// *****************************************************************************************
	//
	void dec() {
		std::unique_lock<std::mutex> lk(mutex);
		--counter;

		if(counter == 0)
			cv.notify_one();
	}

	// *****************************************************************************************
	//
	void dec_notify_all() {
		std::unique_lock<std::mutex> lk(mutex);
		--counter;

		if (counter == 0)
			cv.notify_all();
	}

	// *****************************************************************************************
	//
	void waitForZero() {
		std::unique_lock<std::mutex> lk(mutex);
		cv.wait(lk, [this] {return counter == 0; });
	}
};
#endif

// ************************************************************************************
template<typename T> class TaskManager
{
	using queue_t = vector<T>;

	queue_t q;
	atomic<int64_t> n_tasks;
	atomic<int64_t> n_waiting;
	atomic<int64_t> stage_id;
	atomic_flag is_empty;
	atomic_flag is_completed;

public:
	TaskManager() :
		n_tasks(0),
		n_waiting(0),
		stage_id(0)
	{
		is_empty.clear();
		is_completed.clear();
	}

	void AddTasks(const vector<T>& tasks)
	{
		is_empty.clear();

		q = tasks;
		n_waiting = (int64_t)q.size();
		n_tasks = (int64_t)q.size();
		stage_id++;

		if (n_tasks == 0)
			is_completed.test_and_set();

		stage_id.notify_all();
	}

	bool Pop(T& task)
	{
		while(true)
		{
			int64_t prev_stage_id = stage_id;
			int64_t task_id = n_tasks.fetch_sub(1) - 1;

			if (task_id < 0)
			{
				stage_id.wait(prev_stage_id);
				if (is_completed.test())
					return false;
			}
			else
			{
				task = q[task_id];
				return true;
			}
		}
	}

	void MarkProcessed()
	{
		if (n_waiting.fetch_sub(1) == 1)
		{
			is_empty.test_and_set();
			is_empty.notify_one();
		}
	}

	void Wait()
	{
		is_empty.wait(false);
	}

	void BusyWait()
	{
		while (!is_empty.test())
		{
#if defined(_WIN32)
			_mm_pause();
#elif defined(__aarch64__)
			std::this_thread::yield();
#else
			__builtin_ia32_pause();
#endif
		}
	}
};

// ************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class RegisteringQueue
{
	typedef queue<T, deque<T>> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	int n_elements;
	int size;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;
	condition_variable cv_queue_full;

public:
	// *****************************************************************************************
	//
	RegisteringQueue(int _n_producers, int size = 0)
	{
		Restart(_n_producers, size);
	};

	// *****************************************************************************************
	//
	~RegisteringQueue()
	{};

	// *****************************************************************************************
	//
	void Restart(int _n_producers, int size = 0)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers  = _n_producers;
		n_elements = 0;
		this->size = size;
	}

	// *****************************************************************************************
	//
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *****************************************************************************************
	//
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *****************************************************************************************
	//
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if(!n_producers)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void Push(T data)
	{
		unique_lock<mutex> lck(mtx);

		if (size > 0) {
			cv_queue_full.wait(lck, [this] {return this->n_elements < size; });
		}

//		bool was_empty = n_elements == 0;
		q.push(data);
		++n_elements;

		cv_queue_empty.notify_one();
//		if(was_empty)
//			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void PushRange(vector<T> &data)
	{
		unique_lock<mutex> lck(mtx);
		
		for(auto &x : data)
			q.push(x);
		n_elements += data.size();

//		if (was_empty)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	bool Pop(T &data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this]{return !this->q.empty() || !this->n_producers;}); 

		if(n_elements == 0)
			return false;

		data = q.front();
		q.pop();
		--n_elements;
//		if(n_elements == 0)
//			cv_queue_empty.notify_all();
//		if(size && n_elements + 1 == size)
//			cv_queue_full.notify_all();
		if(size > 0)
			cv_queue_full.notify_one();

		return true;
	}

	// *****************************************************************************************
	//
	uint32_t GetSize()
	{
		unique_lock<mutex> lck(mtx);
		return n_elements;
	}
};




// ************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> 
class SynchronizedPriorityQueue
{
	std::priority_queue<
		std::pair<int, T>, 
		std::vector<std::pair<int, T>>, 
		std::function<bool(const std::pair<int, T>&, const std::pair<int, T>&)>> q;
	
	bool is_completed;
	int n_producers;
	int n_elements;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	// *****************************************************************************************
	//
	SynchronizedPriorityQueue(int _n_producers) : q([](const std::pair<int, T>& a, const std::pair<int, T>& b)->bool { return a.first > b.first; })
	{
		Restart(_n_producers);
	};

	// *****************************************************************************************
	//
	~SynchronizedPriorityQueue()
	{};

	// *****************************************************************************************
	//
	void Restart(int _n_producers)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		n_elements = 0;
	}

	// *****************************************************************************************
	//
	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	// *****************************************************************************************
	//
	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	// *****************************************************************************************
	//
	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void Push(int priority, T data)
	{
		unique_lock<mutex> lck(mtx);
	//	bool was_empty = n_elements == 0;
		q.push(std::make_pair(priority, data));
		++n_elements;

	//	if (was_empty)
		cv_queue_empty.notify_all();
	}

	
	// *****************************************************************************************
	//
	bool Pop(int priority, T &data)
	{
		unique_lock<mutex> lck(mtx);
		cv_queue_empty.wait(lck, [this, priority] {return (!this->q.empty() && q.top().first == priority)  || !this->n_producers; });

		if (n_elements == 0 || q.top().first != priority)
			return false;

		data = q.top().second;	
		q.pop();
		--n_elements;
	//	if (n_elements == 0)
		cv_queue_empty.notify_all();

		return true;
	}

	// *****************************************************************************************
	//
	uint32_t GetSize()
	{
		unique_lock<mutex> lck(mtx);
		return n_elements;
	}
};

#if 0
// ************************************************************************************
// Multithreading queue assuming some special use scheme
//   * Init() can be called when no other thread operates on the queue
template<typename T> class TwoPassQueue
{
	vector<T> q;
	atomic<int> n_elements;
	atomic<int> n_producers;
	atomic<int64_t> id_stage;

public:
	// *****************************************************************************************
	//
	TwoPassQueue(int _n_producers)
	{
		Restart(_n_producers);
	};

	// *****************************************************************************************
	//
	~TwoPassQueue()
	{};

	// *****************************************************************************************
	//
	void Restart(int _n_producers)
	{
		n_producers = _n_producers;
		n_elements = 0;
		id_stage = 0;
	}

	// *****************************************************************************************
	//
	bool IsEmpty()
	{
		return n_elements == 0;
	}

	// *****************************************************************************************
	//
/*	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}*/

	// *****************************************************************************************
	//
	void MarkCompleted()
	{
		int cur_n_producers = n_producers.fetch_sub(1) - 1;

		if (cur_n_producers == 0)
		{
			++id_stage;
			id_stage.notify_all();
		}
	}

	// *****************************************************************************************
	//
//	void Init(vector<T>& data)
	void PushRange(vector<T>& data)
	{
		q = data;
		n_elements = (int) q.size();
		++id_stage;

		id_stage.notify_all();
	}

	// *****************************************************************************************
	//
	bool Pop(T& data)
	{
		int id;

		do
		{
			id = n_elements.fetch_sub(1) - 1;

			if (id < 0)
			{
				if (n_producers == 0)
					return false;

				int64_t old_id_stage = id_stage;
				id_stage.wait(old_id_stage);
			}
		} while (id < 0);

		data = q[id];

		return true;
	}

	// *****************************************************************************************
	//
	uint32_t GetSize()
	{
		return n_elements;
	}
};

#endif
