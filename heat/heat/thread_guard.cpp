#include "thread_guard.h"
//#include "stdafx.h"

thread_guard::thread_guard() {}

void thread_guard::join()
{
    for(int i=0;i < t.size() ;i++){
        if (t[i].joinable()) t[i].join();
    }
}

thread_guard::~thread_guard()
{
    join();
}

int thread_guard::get_available() {
    return (thread::hardware_concurrency()- t.size());
}

bool thread_guard::push(thread in)
{
    _m.lock();
    if (get_available() != 0)
    {
        this->t.push_back(move(in));
        _m.unlock();
        return true;
    }
    else
    {
        _m.unlock();
        return false;
    }
}
