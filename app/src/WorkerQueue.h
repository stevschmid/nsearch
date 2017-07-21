#pragma once

#include <deque>
#include <queue>
#include <thread>

template <class QueueItem>
class WorkerQueue {
public:
  using OnProcessedCallback = std::function< void ( size_t, size_t )  >;

  WorkerQueue( int numWorkers = 1 );
  ~WorkerQueue();

  void Enqueue( QueueItem &queueItem );
  bool Done() const;
  void WaitTillDone();

  void OnProcessed( const OnProcessedCallback &callback ) {
    mProcessedCallbacks.push_back( callback );
  }

protected:
  virtual void Process( const QueueItem &queueItem ) = 0;

private:
  friend class Worker;
  std::deque< std::thread > mWorkers;

  std::condition_variable mCondition;
  std::mutex mQueueMutex;
  std::atomic< bool > mStop;
  std::atomic< int > mWorkingCount;

  std::queue< QueueItem > mQueue;

  size_t mTotalEnqueued;
  size_t mTotalProcessed;
  std::deque< OnProcessedCallback > mProcessedCallbacks;

  void WorkerLoop();
};

template <class QueueItem>
WorkerQueue< QueueItem >::WorkerQueue( int numWorkers )
  : mStop( false ), mWorkingCount( 0 ), mTotalEnqueued( 0 ), mTotalProcessed( 0 )
{
  numWorkers = numWorkers <= 0 ? std::thread::hardware_concurrency() : numWorkers;

  for( int i = 0; i < numWorkers; i++ ) {
    mWorkers.push_back( std::thread( [this] { this->WorkerLoop(); } ) );
  }
}

template <class QueueItem>
WorkerQueue< QueueItem > ::~WorkerQueue() {
  mStop = true;
  mCondition.notify_all();
  for( auto &worker : mWorkers ) {
    if( worker.joinable() ) {
      worker.join();
    }
  }
}

template <class QueueItem>
bool WorkerQueue< QueueItem >::Done() const {
  return mWorkingCount == 0 && mQueue.empty();
}

template <class QueueItem>
void WorkerQueue< QueueItem >::WaitTillDone() {
  while( !Done() ) {
    std::this_thread::sleep_for( std::chrono::milliseconds( 50 ) );
  }
}

template<class QueueItem>
void WorkerQueue< QueueItem >::Enqueue( QueueItem &queueItem )
{
  {
    std::unique_lock< std::mutex > lock( mQueueMutex );
    mTotalEnqueued++;
    mQueue.push( std::move( queueItem ) );
  }

  // Send one worker to work
  mCondition.notify_one();
}

template<class QueueItem>
void WorkerQueue< QueueItem >::WorkerLoop() {
  QueueItem queueItem;

  while( true )
  {
    { // acquire lock
      std::unique_lock< std::mutex > lock( mQueueMutex );

      while( !mStop && mQueue.empty() )
        mCondition.wait( lock );

      if( mStop )
        break;

      queueItem = std::move( mQueue.front() );
      mQueue.pop();

      mWorkingCount++;
    } // release lock

    Process( queueItem );

    { // acquire lock
      std::unique_lock< std::mutex > lock( mQueueMutex );
      mTotalProcessed++;
      mWorkingCount--;

      for( auto &cb : mProcessedCallbacks ) {
        cb( mTotalProcessed, mTotalEnqueued );
      }
    } // release lock
  }
}
