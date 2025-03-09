// timeUtils.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "TimeUtils.h"

std::string TimeUtils::timeStampToString(std::time_t inTS, const char* sFormatString)
{
    assert(s_bTested);
    std::tm tmpTM = timeStampToTM(inTS);
    char buffer[80];
    std::strftime(buffer, sizeof(buffer), sFormatString, &tmpTM);
    return buffer;
}
