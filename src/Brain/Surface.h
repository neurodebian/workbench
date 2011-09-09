
#ifndef __SURFACE_H__
#define __SURFACE_H__

/*LICENSE_START*/ 
/* 
 *  Copyright 1995-2002 Washington University School of Medicine 
 * 
 *  http://brainmap.wustl.edu 
 * 
 *  This file is part of CARET. 
 * 
 *  CARET is free software; you can redistribute it and/or modify 
 *  it under the terms of the GNU General Public License as published by 
 *  the Free Software Foundation; either version 2 of the License, or 
 *  (at your option) any later version. 
 * 
 *  CARET is distributed in the hope that it will be useful, 
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of 
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details. 
 * 
 *  You should have received a copy of the GNU General Public License 
 *  along with CARET; if not, write to the Free Software 
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 
 * 
 */ 

#include <vector>

#include "SurfaceFile.h"

namespace caret {
    
    class BoundingBox;
    class Brain;
    class BrainStructure;
    class ModelDisplayController;
    class ModelDisplayControllerSurface;

    /**
     * Maintains view of some type of object.
     */
    class Surface : public SurfaceFile {
        
    public:
        Surface();
        
        ~Surface();

        Surface(const Surface& s);
        
        Surface& operator=(const Surface& s);
        
        void setBrainStructure(BrainStructure* brainStructure);
        
        ModelDisplayController* getModelController();
        
        AString getNameForGUI(bool includeStructureFlag) const;
        
        Brain* getBrain();
        
        BrainStructure* getBrainStructure();
        
        const BrainStructure* getBrainStructure() const;
        
        void getBounds(BoundingBox& boundingBoxOut) const;
        
    private:
        void initializeMemberSurface();
        
        void copyHelperSurface(const Surface& s);
        
        BrainStructure* brainStructure;
        
        ModelDisplayControllerSurface* surfaceController;
        
        bool defaultScalingInitializedFlag;
        
        std::vector<float> normalVectors;
    };

} // namespace

#endif // __SURFACE_H__
