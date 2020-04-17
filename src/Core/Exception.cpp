#include <Endas/Core/Exception.hpp>
#include <sstream>
 
using namespace std;
using namespace endas;


template <class MsgFormat>
inline string formatFileLineMsg(const char* file, int line, MsgFormat msg)
{
    stringstream ss;
    ss << file << ":" << line << ": ";
    msg(ss);
    return ss.str();
}


NotSupportedError::NotSupportedError(const char* msg, const char* file, int line)
: std::logic_error(formatFileLineMsg(file, line, [=](ostream& _msg) { _msg << msg; }))
{ }


NotImplementedError::NotImplementedError(const char* file, int line)
: std::logic_error(formatFileLineMsg(file, line, [](ostream& msg) { msg << "Not implemented"; }))
{ }


AssertionError::AssertionError(const char* expr, const char* file, int line)
: std::logic_error(formatFileLineMsg(file, line, [=](ostream& msg) { msg << "Assertion '" << expr << "' failed"; }))
{ }
