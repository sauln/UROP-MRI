


#ifndef DEBUG_STREAM_H
#define DEBUG_STREAM_H
 #include <iostream>
 #include <string>
 class DebugStream
 {
	     std::ostream& m_os;
     bool m_debug;
	
		 public:
	     DebugStream(std::ostream& os) :
				         m_os(os),
				         m_debug(true)
				     { }
			
				     bool debug() const { return m_debug; }
			   void setDebug(bool b) { m_debug = b; }
			
				     std::ostream& os() { return m_os; }
			
				     DebugStream& operator<< (bool val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			    DebugStream& operator<< (short val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
		    DebugStream& operator<< (unsigned short val)
			     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			    DebugStream& operator<< (int val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			     DebugStream& operator<< (unsigned int val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			    DebugStream& operator<< (long val)
				     {
			         if (m_debug) m_os << val;
				         return *this;
				     }
			     DebugStream& operator<< (unsigned long val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			    DebugStream& operator<< (float val)
				     {
				         if (m_debug) m_os << val;
				        return *this;
				     }
			    DebugStream& operator<< (double val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			     DebugStream& operator<< (long double val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			     DebugStream& operator<< (const void* val)
				     {
				         if (m_debug) m_os << val;
				         return *this;
				     }
			
				     DebugStream& operator<< (std::streambuf* sb)
				     {
				         if (m_debug) m_os << sb;
				         return *this;
				     }
		
				     DebugStream& operator<< (std::ostream& (*pf)(std::ostream&))
				     {
				         if (m_debug) m_os << pf;
				         return *this;
				     }
			     DebugStream& operator<< (std::ios& (*pf)(std::ios&))
				     {
				        if (m_debug) m_os << pf;
				         return *this;
				     }
			     DebugStream& operator<< (std::ios_base& (*pf)(std::ios_base&))
				     {
				         if (m_debug) m_os << pf;
				         return *this;
				     }
			
				     DebugStream& put(char c)
				   {
				        if (m_debug) m_os.put(c);
				       return *this;
				     }
			
				     DebugStream& write(const char* s, std::streamsize n)
				     {
				         if (m_debug) m_os.write(s, n);
				         return *this;
				     }
			
				    friend DebugStream& operator<< (DebugStream& out, char c);
			     friend DebugStream& operator<< (DebugStream& out, signed char c);
			    friend DebugStream& operator<< (DebugStream& out, unsigned char c);
			
			    friend DebugStream& operator<< (DebugStream& out, const char* c);
			     friend DebugStream& operator<< (DebugStream& out, const signed char* c);
			     friend DebugStream& operator<< (DebugStream& out, const unsigned char* c);
			 };

 template < class T >
 DebugStream& operator<< (DebugStream& out, const T& t)
 {
	    if (out.debug()) out.os() << t;
	    return out;
	 }

 #endif // DEBUG_STREAM_H