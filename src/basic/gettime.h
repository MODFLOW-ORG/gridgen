
#pragma once
#ifndef _GETTIME_H_
#define _GETTIME_H_

//-----------------------------------------------------------------------------

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include <stdio.h>
#include <stdlib.h>

//------------------------------------------------------------------------
inline double getTime() //in millisecond
{
#ifdef WIN32 || _WIN64
    return timeGetTime();
#else 
    //assuming unix-type systems
    //timezone tz;
    timeval  tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec*1000000+tv.tv_usec)*1.0/1000;
#endif
}

#endif //_GETTIME_H_

