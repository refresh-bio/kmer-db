#ifndef _QUEUE_H
#define _QUEUE_H
// Generic multithreading queues

//#include "../asm_common/defs.h"
#include <queue>
#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;

// ************************************************************************************
// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class CRegisteringQueue
{
	typedef queue<T, list<T>> queue_t;

	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;

public:
	CRegisteringQueue(int _n_producers) 
	{
		Restart(_n_producers);
	};

	~CRegisteringQueue() 
	{};

	void Restart(int _n_producers)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers  = _n_producers;
		n_elements = 0;
	}

	bool IsEmpty() 
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	bool IsCompleted() 
	{
		lock_guard<mutex> lck(mtx);

		return n_elements == 0 && n_producers == 0;
	}

	void MarkCompleted() 
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;

		if(!n_producers)
			cv_queue_empty.notify_all();
	}

	void Push(T data) 
	{
		unique_lock<mutex> lck(mtx);
		bool was_empty = n_elements == 0;
		q.push(data);
		++n_elements;

		if(was_empty)
			cv_queue_empty.notify_all();
	}
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

	uint32_t GetSize() 
	{
		return n_elements;
	}
};

#endif

// EOF
