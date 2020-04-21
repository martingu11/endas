#include <Endas/Core/Profiling.hpp>

#include <chrono>
#include <algorithm>

#include <map>
#include <unordered_map>
#include <list>
#include <cstring>
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace endas;

#if !ENDAS_PROFILING_DISABLED

struct PerfScope
{
    const char* key;
    perfclock_t::time_point start;
    perfclock_t::duration duration;
    list<PerfScope> children;
    PerfScope* parent;

    map<const char*, perfclock_t::duration> records; 

    PerfScope() : PerfScope("MAIN", nullptr) { }

    PerfScope(const char* k, PerfScope* par = nullptr)
    : key(k), start(perfclock_t::now()), duration(0), parent(par)
    { }
};

static PerfScope globRootScope;
static PerfScope* globTopScope;


inline PerfScope& topScope()
{
    if (globTopScope == nullptr) globTopScope = &globRootScope;
    return *globTopScope;
}


endas::detail::NewPerfScope::NewPerfScope(const char* key)
{
    auto& top = topScope();
    
    auto childit = std::find_if(top.children.begin(), top.children.end(), 
        [=](const PerfScope& s) { return strcmp(s.key, key) == 0; });

    if (childit == top.children.end())
    {
        topScope().children.emplace_back(key, &topScope());
        globTopScope = &topScope().children.back();
    }
    else
    {
        globTopScope = &(*childit);
        globTopScope->start = perfclock_t::now();
    }
}


endas::detail::NewPerfScope::~NewPerfScope()
{
    topScope().duration+= perfclock_t::now() - topScope().start;
    globTopScope = topScope().parent;
}


void endas::detail::recordTime(const char* key, perfclock_t::time_point start, perfclock_t::time_point end)
{
    topScope().records[key]+= (end - start);
}



/*size_t getLongestKey(const PerfScope& ps)
{
    size_t l = max(strlen(ps.key), 
        max_element(ps.children.begin(), ps.children.end(), [](const PerfScope& ps   ) )
  
}*/

inline string indent(int n, const string& s)
{
    string str;
    for (int i = 0; i != n; i++) str+= "  ";
    str+= s;
    return str;
}

inline string formatDuration(perfclock_t::duration d)
{
    int64_t us = chrono::duration_cast<chrono::microseconds>(d).count();

    stringstream os;
    
    if (us < 1000) os << us << "us";
    else if (us < 1000000) os << fixed << setprecision(3) <<  (us / 1.0e3) << "ms";
    else
    {
        double seconds = us / 1.0e6;
        if (seconds < 60.0) os << fixed << setprecision(3) << seconds << "s";
        else 
        {
            int minutes = (int)(seconds / 60);
            seconds = seconds - (minutes * 60);
            os << minutes << "m " << fixed << setprecision(3) << seconds << "s";
        }
    }

    return os.str();
}

inline void formatRelDuration(ostream& os, perfclock_t::duration d, perfclock_t::duration parent)
{
    if (parent.count() == 0) return;
    double drel = (d * 100.0) / parent;
    double drel2 = (d * 100.0) / globRootScope.duration;

    os << "(";
    os << fixed << setw(5) << right << setprecision(1) << drel << "% / ";
    os << fixed << setw(5) << right << setprecision(1) << drel2 << "%";
    os << ")";
}


void printScope(std::ostream& os, const PerfScope& scope, int level, int col1width, int maxNesting)
{
    // Print scope name and total execution time
    os << setw(col1width) << left << setfill('.') << indent(level, scope.key);
    os << setfill(' ');

    os << ": " << setw(10) << left << formatDuration(scope.duration);
    os << " ";
    if (scope.parent) formatRelDuration(os, scope.duration, scope.parent->duration);
    os << endl;

    if (level+1 <= maxNesting)
    {

        // Then all time records for this scope
        for (auto&& rec : scope.records)
        {
            os << setw(col1width) << left << setfill('.') << indent(level+1, rec.first);
            os << setfill(' ');

            os << ": " << setw(10) << left << formatDuration(rec.second);
            os << " ";
            formatRelDuration(os, rec.second, scope.duration);
            os << endl;
        }

        // And finally nested scopes
        
        for (auto&& child : scope.children)
        {
            printScope(os, child, level+1, col1width, maxNesting);   
        }
    }
    
};


#endif



void endas::profilerClear()
{
#if !ENDAS_PROFILING_DISABLED
    globRootScope = PerfScope();
    globTopScope = nullptr;
#endif
}


void endas::profilingSummary(std::ostream& os, int maxNesting)
{
#if !ENDAS_PROFILING_DISABLED

    globRootScope.duration = perfclock_t::now() - globRootScope.start;

    int col1width = 40;
    printScope(os, topScope(), 0, col1width, maxNesting);
#endif    
}