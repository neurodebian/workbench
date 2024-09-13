#ifndef __CARET_MAPPABLE_DATA_FILE_CLUSTER_FINDER_H__
#define __CARET_MAPPABLE_DATA_FILE_CLUSTER_FINDER_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2024 Washington University School of Medicine
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*LICENSE_END*/


#include <cstdint>
#include <memory>
#include <vector>

#include "BrainordinateCluster.h"
#include "CaretObject.h"
#include "CaretResult.h"
#include "Vector3D.h"

namespace caret {

    class CaretMappableDataFile;
    class VolumeFile;
    
    class CaretMappableDataFileClusterFinder : public CaretObject {
//        class Cluster {
//        public:
//            Cluster(const AString  roiName,
//                         const int32_t  labelIndex,
//                         const Vector3D roiCenterXYZ,
//                         const int64_t  numberOfBrainordinatesInROI)
//            : m_roiName(roiName),
//            m_labelIndex(labelIndex),
//            m_roiCenterXYZ(roiCenterXYZ),
//            m_numberOfBrainordinatesInROI(numberOfBrainordinatesInROI)
//            { }
//            
//            const AString m_roiName;
//            const int32_t m_labelIndex;
//            const Vector3D m_roiCenterXYZ;
//            const int64_t m_numberOfBrainordinatesInROI;
//        };
            
    public:
        /**
         * Mode for searchinhg
         */
        enum class FindMode {
            VOLUME_LABEL
        };
        
        CaretMappableDataFileClusterFinder(const FindMode findMode,
                                           const CaretMappableDataFile* mapFile,
                                           const int32_t mapIndex);
        
        virtual ~CaretMappableDataFileClusterFinder();
        
        CaretMappableDataFileClusterFinder(const CaretMappableDataFileClusterFinder&) = delete;

        CaretMappableDataFileClusterFinder& operator=(const CaretMappableDataFileClusterFinder&) = delete;
        
        std::unique_ptr<CaretResult> findClusters();
        
        const std::vector<BrainordinateCluster>& getClusters() const;
        
        AString getClustersInFormattedString() const;
        
        // ADD_NEW_METHODS_HERE

    private:
        std::unique_ptr<CaretResult> findLabelVolumeClusters(const VolumeFile* volumeFile);
        
        const FindMode m_findMode;
        
        const CaretMappableDataFile* m_mapFile;
        
        const int32_t m_mapIndex;
        
        std::vector<BrainordinateCluster> m_clusters;
        
//        float m_valueLow = 0.0;
//        
//        float m_valueHigh = 0.0;
        
//        std::vector<int8_t> m_voxelSearchedFlag;
//        
//        VolumeSpace m_volumeSpace;
        
        // ADD_NEW_MEMBERS_HERE

    };
    
#ifdef __CARET_MAPPABLE_DATA_FILE_CLUSTER_FINDER_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __CARET_MAPPABLE_DATA_FILE_CLUSTER_FINDER_DECLARE__

} // namespace
#endif  //__CARET_MAPPABLE_DATA_FILE_CLUSTER_FINDER_H__
