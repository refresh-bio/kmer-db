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

using namespace std;



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

		cv.notify_one();
	}

	// *****************************************************************************************
	//
	void dec_notify_all() {
		std::unique_lock<std::mutex> lk(mutex);
		--counter;

		cv.notify_all();
	}

	// *****************************************************************************************
	//
	void waitForZero() {
		std::unique_lock<std::mutex> lk(mutex);
		cv.wait(lk, [this] {return counter == 0; });
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
	uint32_t n_elements;
	int size;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

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
			cv_queue_empty.wait(lck, [this] {return this->n_elements < size; });
		}

		bool was_empty = n_elements == 0;
		q.push(data);
		++n_elements;

		if(was_empty)
			cv_queue_empty.notify_all();
	}

	// *****************************************************************************************
	//
	void PushRange(vector<T> &data)
	{
		unique_lock<mutex> lck(mtx);
		bool was_empty = n_elements == 0;
		
		for(auto &x : data)
			q.push(x);
		n_elements += data.size();

		if (was_empty)
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
		if(n_elements == 0)
			cv_queue_empty.notify_all();

		return true;
	}

	// *****************************************************************************************
	//
	uint32_t GetSize()
	{
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
	uint32_t n_elements;

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

		if (n_elements == 0)
			return false;

		auto entry = q.top();
		data = entry.second;
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
		return n_elements;
	}
};