#pragma once

#include <deque>
#include <queue>
#include <thread>

class Worker;

class ThreadPool {
private:
  friend class Worker;
  std::deque< std::thread > mWorkers;

  std::condition_variable mCondition;
  std::mutex mQueueMutex;
  std::atomic< bool > mStop;
  std::atomic< int > mWorkingCount;

  std::queue< std::function< void() > > mQueue;

public:
  ThreadPool( int numWorkers = std::thread::hardware_concurrency() );
  ~ThreadPool();

  template<class F>
  void Enqueue(F f);

  bool Done();
};

class Worker {
private:
  ThreadPool &mThreadPool;

public:
  Worker( ThreadPool &threadPool );
  void operator()();
};

ThreadPool::ThreadPool( int numWorkers )
  : mStop( false ), mWorkingCount( 0 )
{
  for( int i = 0; i < numWorkers; i++ ) {
    mWorkers.push_back( std::thread( Worker( *this ) ) );
  }
}

ThreadPool::~ThreadPool() {
  mStop = true;
  mCondition.notify_all();
  for( auto &worker : mWorkers ) {
    if( worker.joinable() ) {
      worker.join();
    }
  }
}

bool ThreadPool::Done() {
  std::lock_guard< std::mutex > lock( mQueueMutex );
  return mWorkingCount == 0 && mQueue.empty();
}

template<class F>
void ThreadPool::Enqueue(F f)
{
  {
    std::unique_lock< std::mutex > lock( mQueueMutex );
    mQueue.push( std::function< void() >( f ) );
  }

  // Send one worker to work
  mCondition.notify_one();
}

Worker::Worker( ThreadPool &threadPool )
  : mThreadPool( threadPool )
{
}

void Worker::operator()() {
  std::function< void() > task;

  while( true )
  {
    { // acquire lock
      std::unique_lock< std::mutex > lock( mThreadPool.mQueueMutex );

      while( !mThreadPool.mStop && mThreadPool.mQueue.empty() )
        mThreadPool.mCondition.wait( lock );

      if( mThreadPool.mStop )
        break;

      task = mThreadPool.mQueue.front();
      mThreadPool.mQueue.pop();

      mThreadPool.mWorkingCount++;
    } // release lock

    task();

    { // acquire lock
      std::unique_lock< std::mutex > lock( mThreadPool.mQueueMutex );
      mThreadPool.mWorkingCount--;
    }
  }
}
