/*
 *  concurrent_queue.h
 *  rda
 *
 *  Created by Alexander K. Hudek on 2008-08-05.
 *  Copyright 2008 University of Waterloo. All rights reserved.
 *
 */
#ifndef CONCURRENT_QUEUE__
#define CONCURRENT_QUEUE__
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <queue>

template< typename DT >
class concurrent_queue {
protected:
	std::queue<DT> my_queue;
	mutable boost::mutex my_mutex;
	mutable boost::condition cond;
	bool finished, finish_when_empty;
	
public:
	typedef DT const & const_reference;
	
	concurrent_queue() : finished(false), finish_when_empty(false) {}
	
	void push_back( DT const &v ) {
		boost::mutex::scoped_lock lock(my_mutex);
		bool was_empty = my_queue.empty();
		my_queue.push(v);
		lock.unlock();
		if( was_empty ) cond.notify_one();
	}
	
	bool empty() const {
		boost::mutex::scoped_lock lock(my_mutex);
		return my_queue.empty();
	}
	
	void shutdown_when_empty() {
		boost::mutex::scoped_lock lock(my_mutex);
		if( my_queue.empty() ) {
			lock.unlock();
			shutdown();
		} else finish_when_empty = true;
	}

	// shuts the queue down and all threads get false instead of data
	void shutdown() {
		boost::mutex::scoped_lock lock(my_mutex);
		finished = true;
		lock.unlock();
		cond.notify_all();		
	}
	
	// Gets data off the queue. This will wait until data is available.
	// If it returns false, then release has been called.
	bool get_data( DT &r ) {
		boost::mutex::scoped_lock lock(my_mutex);

		// wait for data or for the release notification
		while( !finished && my_queue.empty() ) {
			cond.wait(lock);
		}
		
		if( finished ) return false;
		
		r = my_queue.front();
		my_queue.pop();
		
		// finish up if we are done
		if( my_queue.empty() && finish_when_empty ) {
			lock.unlock();
			shutdown();
		}
		return true;
	}
	
};

#endif
