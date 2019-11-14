#ifndef THREAD_GUARD_H
#define THREAD_GUARD_H

#include <thread>
#include <vector>
#include <mutex>

using namespace std;

class thread_guard
{
private:
    vector<thread> t;
    mutex _m;
    //int get_available();
public:
    thread_guard();
    bool push(thread t);
    void join();
    ~thread_guard();
    thread_guard (thread_guard& in) = delete;
    thread_guard& operator = (thread_guard& in) = delete;
	int get_available();
};


#endif //THREAD_GUARD_H
