#ifndef __TILE_TABS_LAYOUT_MANUAL_CONFIGURATION_H__
#define __TILE_TABS_LAYOUT_MANUAL_CONFIGURATION_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2019 Washington University School of Medicine
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



#include <memory>

#include "TileTabsLayoutBaseConfiguration.h"
#include "TileTabsLayoutConfigurationTypeEnum.h"

namespace caret {

    class BrowserTabContent;
    class TileTabsLayoutGridConfiguration;
    class TileTabsBrowserTabGeometry;
    
    class TileTabsLayoutManualConfiguration : public TileTabsLayoutBaseConfiguration {
        
    public:
        static TileTabsLayoutManualConfiguration* newInstanceFromGridLayout(TileTabsLayoutGridConfiguration* gridLayout,
                                                                            const TileTabsLayoutConfigurationTypeEnum::Enum gridMode,
                                                                            const std::vector<int32_t>& tabIndices);
        
        TileTabsLayoutManualConfiguration();
        
        virtual ~TileTabsLayoutManualConfiguration();
        
        TileTabsLayoutManualConfiguration(const TileTabsLayoutManualConfiguration& obj);

        TileTabsLayoutManualConfiguration& operator=(const TileTabsLayoutManualConfiguration& obj);
        
        virtual void copy(const TileTabsLayoutBaseConfiguration& rhs) override;
        
        virtual TileTabsLayoutBaseConfiguration* newCopyWithNewUniqueIdentifier() const override;
        
        virtual AString toString() const override;
        
        void addTabInfo(TileTabsBrowserTabGeometry* tabInfo);
        
        virtual int32_t getNumberOfTabs() const override;
        
        TileTabsBrowserTabGeometry* getTabInfo(const int32_t index);
        
        const TileTabsBrowserTabGeometry* getTabInfo(const int32_t index) const;
        
        virtual TileTabsLayoutManualConfiguration* castToManualConfiguration() override;
        
        virtual const TileTabsLayoutManualConfiguration* castToManualConfiguration() const override;
        
        int32_t getWindowAnnotationsStackingOrder() const;
        
        void setWindowAnnotationsStackingOrder(const int32_t stackingOrder);
        
        void initializeManualConfiguration();
        
        // ADD_NEW_METHODS_HERE

    protected:
        virtual void decodeFromXMLString(QXmlStreamReader& xml,
                                         const AString& rootElementText) override;
        
        virtual void encodeInXMLString(AString& xmlTextOut) const override;
        
    private:
        void correctOldStackingOrder();
        
        void copyHelperTileTabsLayoutManualConfiguration(const TileTabsLayoutManualConfiguration& obj);

        std::vector<std::unique_ptr<TileTabsBrowserTabGeometry>> m_tabInfo;
        
        int32_t m_windowAnnotationsStackingOrder = -1000;
        

        // ADD_NEW_MEMBERS_HERE

        
        static const AString s_rootElementName;
        static const AString s_rootElementAttributeVersion;
        static const AString s_rootElementAttributeName;
        static const AString s_rootElementAttributeUniqueID;

        static const AString s_rootElementAttributeValueVersionOne;
        
        static const AString s_tabInfoElementName;
        static const AString s_tabInfoAttributeDisplayStatus;
        static const AString s_tabInfoAttributeTabIndex;
        static const AString s_tabInfoAttributeMinX;
        static const AString s_tabInfoAttributeMaxX;
        static const AString s_tabInfoAttributeMinY;
        static const AString s_tabInfoAttributeMaxY;
        static const AString s_tabInfoAttributeStackingOrder;
        static const AString s_tabInfoAttributeBackground;

        static const AString s_windowInfoElementName;
        static const AString s_windowInfoAttributeStackingOrder;
        
        friend class TileTabsLayoutBaseConfiguration;
    };
    
#ifdef __TILE_TABS_LAYOUT_MANUAL_CONFIGURATION_DECLARE__
    const AString TileTabsLayoutManualConfiguration::s_rootElementName = "TileTabsManualLayout";
    const AString TileTabsLayoutManualConfiguration::s_rootElementAttributeVersion = "Version";
    const AString TileTabsLayoutManualConfiguration::s_rootElementAttributeName = "Name";
    const AString TileTabsLayoutManualConfiguration::s_rootElementAttributeUniqueID = "UniqueID";
    
    const AString TileTabsLayoutManualConfiguration::s_rootElementAttributeValueVersionOne = "1";

    const AString TileTabsLayoutManualConfiguration::s_tabInfoElementName = "TabInfo";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeDisplayStatus = "DisplayStatus";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeTabIndex = "TabIndex";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeMinX = "MinX";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeMaxX = "MaxX";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeMinY = "MinY";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeMaxY = "MaxY";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeStackingOrder = "StackingOrder";
    const AString TileTabsLayoutManualConfiguration::s_tabInfoAttributeBackground = "Background";

    const AString TileTabsLayoutManualConfiguration::s_windowInfoElementName            = "WindowInfo";
    const AString TileTabsLayoutManualConfiguration::s_windowInfoAttributeStackingOrder = "StackingOrder";
    
#endif // __TILE_TABS_LAYOUT_MANUAL_CONFIGURATION_DECLARE__

} // namespace
#endif  //__TILE_TABS_LAYOUT_MANUAL_CONFIGURATION_H__
