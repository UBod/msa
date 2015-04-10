#ifndef RCLUSTALWMAIN_H_
#define RCLUSTALWMAIN_H_

typedef std::vector<std::string> UserArgs;

#include "general/clustalw.h"
#include <Rcpp.h>

using namespace std;

namespace clustalw {
	class RClustalWMain {
		public:
			RClustalWMain();
			~RClustalWMain();
			void run(UserArgs args, ClustalWInput *input, ClustalWOutput *output);
	};
}
#endif
