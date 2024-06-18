#ifndef __BACKGROUND_AND_FOREGROUND_COLORS_H__
#define __BACKGROUND_AND_FOREGROUND_COLORS_H__

/*LICENSE_START*/
/*
 *  Copyright (C) 2016 Washington University School of Medicine
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

#include <array>
#include "CaretObject.h"



namespace caret {

    class BackgroundAndForegroundColors : public CaretObject {
        
    public:
        static std::array<float, 3> toFloatRGB(const uint8_t byteRGB[3]);
        
        BackgroundAndForegroundColors();
        
        virtual ~BackgroundAndForegroundColors();
        
        BackgroundAndForegroundColors(const BackgroundAndForegroundColors& obj);

        BackgroundAndForegroundColors& operator=(const BackgroundAndForegroundColors& obj);

        bool operator==(const BackgroundAndForegroundColors& obj) const;
        
        void reset();

        // ADD_NEW_METHODS_HERE

        void getColorBackgroundWindow(uint8_t colorBackground[3]) const;
        
        void setColorBackgroundWindow(const uint8_t colorBackground[3]);
        
        void getColorForegroundWindow(uint8_t colorForeground[3]) const;
        
        void setColorForegroundWindow(const uint8_t colorForeground[3]);
        

        void getColorForegroundAllView(uint8_t colorForeground[3]) const;
        
        void setColorForegroundAllView(const uint8_t colorForeground[3]);
        
        void getColorBackgroundAllView(uint8_t colorForeground[3]) const;
        
        void setColorBackgroundAllView(const uint8_t colorForeground[3]);
        
        
        void getColorForegroundChartView(uint8_t colorForeground[3]) const;
        
        void setColorForegroundChartView(const uint8_t colorForeground[3]);
        
        void getColorBackgroundChartView(uint8_t colorForeground[3]) const;
        
        void setColorBackgroundChartView(const uint8_t colorForeground[3]);
        
        
        void getColorForegroundSurfaceView(uint8_t colorForeground[3]) const;
        
        void setColorForegroundSurfaceView(const uint8_t colorForeground[3]);
        
        void getColorBackgroundSurfaceView(uint8_t colorForeground[3]) const;
        
        void setColorBackgroundSurfaceView(const uint8_t colorForeground[3]);
        
        
        void getColorForegroundVolumeView(uint8_t colorForeground[3]) const;
        
        void setColorForegroundVolumeView(const uint8_t colorForeground[3]);
        
        void getColorBackgroundVolumeView(uint8_t colorForeground[3]) const;
        
        void setColorBackgroundVolumeView(const uint8_t colorForeground[3]);
        
        void getColorChartMatrixGridLines(uint8_t colorChartMatrixGridLines[3]) const;
        
        void setColorChartMatrixGridLines(const uint8_t colorChartMatrixGridLines[3]);
        
        void getColorChartHistogramThreshold(uint8_t colorChartMatrixGridLines[3]) const;
        
        void setColorChartHistogramThreshold(const uint8_t colorChartMatrixGridLines[3]);
        
        void getColorForegroundHistologyView(uint8_t colorForeground[3]) const;
        
        void setColorForegroundHistologyView(const uint8_t colorForeground[3]);
        
        void getColorBackgroundHistologyView(uint8_t colorForeground[3]) const;
        
        void setColorBackgroundHistologyView(const uint8_t colorForeground[3]);
        
        void getColorForegroundMediaView(uint8_t colorForeground[3]) const;
        
        void setColorForegroundMediaView(const uint8_t colorForeground[3]);
        
        void getColorBackgroundMediaView(uint8_t colorForeground[3]) const;
        
        void setColorBackgroundMediaView(const uint8_t colorForeground[3]);
        
    private:
        void copyHelperBackgroundAndForegroundColors(const BackgroundAndForegroundColors& obj);

        void setColor(uint8_t color[3],
                      const uint8_t red,
                      const uint8_t green,
                      const uint8_t blue);
        
        uint8_t m_colorBackgroundWindow[3];
        
        uint8_t m_colorForegroundWindow[3];
        
        uint8_t m_colorForegroundAll[3];
        
        uint8_t m_colorBackgroundAll[3];
        
        uint8_t m_colorForegroundChart[3];
        
        uint8_t m_colorBackgroundChart[3];
        
        uint8_t m_colorForegroundSurface[3];
        
        uint8_t m_colorBackgroundSurface[3];
        
        uint8_t m_colorForegroundVolume[3];
        
        uint8_t m_colorBackgroundVolume[3];
        
        uint8_t m_colorChartMatrixGridLines[3];
        
        uint8_t m_colorChartHistogramThreshold[3];
        
        uint8_t m_colorForegroundHistology[3];
        
        uint8_t m_colorBackgroundHistology[3];
        
        uint8_t m_colorForegroundMedia[3];
        
        uint8_t m_colorBackgroundMedia[3];
        
        // ADD_NEW_MEMBERS_HERE

        friend class BackgroundAndForegroundColorsSceneHelper;
        friend class CaretPreferences;
    };
    
#ifdef __BACKGROUND_AND_FOREGROUND_COLORS_DECLARE__
    // <PLACE DECLARATIONS OF STATIC MEMBERS HERE>
#endif // __BACKGROUND_AND_FOREGROUND_COLORS_DECLARE__

} // namespace
#endif  //__BACKGROUND_AND_FOREGROUND_COLORS_H__
