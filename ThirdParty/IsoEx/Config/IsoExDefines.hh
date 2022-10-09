#ifndef ISOEXDLLEXPORT
	#ifdef WIN32
		#ifdef ISOEXDLL
			#ifdef USEISOEX
				#define ISOEXDLLEXPORT __declspec(dllimport)
				#define ISOEXDLLEXPORTONLY 
			#else
				#define ISOEXDLLEXPORT __declspec(dllexport)
				#define ISOEXDLLEXPORTONLY __declspec(dllexport)
			#endif
		#else		
			#define ISOEXDLLEXPORT	
			#define ISOEXDLLEXPORTONLY
		#endif
	#else
		#define ISOEXDLLEXPORT
		#define ISOEXDLLEXPORTONLY
	#endif
#endif
