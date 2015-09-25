#ifndef __WRITEUTIL__
#define __WRITEUTIL__

#include <iostream>

class WriteUtil
{
 public:

  static void OpenOutputDst(const std::string& dst_name);

  static void CloseOutputDst();

 private:

  WriteUtil();
  ~WriteUtil();
  WriteUtil(WriteUtil const&);
  WriteUtil& operator=(WriteUtil const&);

};

#endif
