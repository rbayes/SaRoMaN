// ----------------------------------------------------------------------------
///  \file     MindMaterialsList.h
///  \brief    Definition of materials used in geometry 
///
///  \author   J Martin-Albo <jmalbos@ific.uv.es>
///  \date     4 Feb, 2009
///  \version  $Id: MindMaterialsList.h 263 2009-06-03 16:49:36Z jmalbos $
///
///  Copyright (c) 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __MATERIALS_LIST__
#define __MATERIALS_LIST__


/// Definition of the materials used in the detector geometry. The class 
/// contains only a static method that must be invoked before the detector 
/// construction.

class MindMaterialsList
{
public:
  /// This static method must be invoked before the detector construction.
  static void DefineMaterials();
  
private:
  /// Constructor (hidden)
  MindMaterialsList() {}
  /// Destructor (hidden)
  ~MindMaterialsList() {}
};

#endif
