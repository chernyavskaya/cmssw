#ifndef EventCategory_H
#define EventCategory_H

/** \class EventCategory
 *  No description available.
 *
 *  $Date: 2011/08/05 08:54:56 $
 *  $Revision: 1.1 $
 *  \author L. Quertenmont P. Silva
 */



#include <utility>
#include <vector>
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuSummaryHandler.h"
#include "CMGTools/HtoZZ2l2nu/interface/ZZ2l2nuPhysicsEvent.h"

const TString ZZ2l2nuCategoryLabel[] = {"eq0jets","eq1jets","geq2jets", "vbf"};

class EventCategory {
public:
  enum CategoryType { EQ0JETS, EQ1JETS, GEQ2JETS, VBF };

 
  /// Constructor
  EventCategory();
  /// Destructor
  virtual ~EventCategory();

  int Get(const PhysicsEvent_t& phys);
  TString GetLabel(int CategoryType);
  TString GetLabel(const PhysicsEvent_t& phys);

};
#endif

