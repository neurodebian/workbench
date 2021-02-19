/*LICENSE_START*/
/*
 *  Copyright (C) 2020  Washington University School of Medicine
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

#include "PaletteNew.h"

#include "CaretAssert.h"
#include "CaretException.h"
#include "CaretLogger.h"

#include "ColorFunctions.h"
#include "MathFunctions.h"

#include <algorithm>

using namespace std;
using namespace caret;

namespace
{
    int searchLowSide(const float scalar, const vector<PaletteNew::ScalarColor>& controlPoints, int low, int high)
    {
        while (low < high - 1)
        {
            int guess = (high + low) / 2;//simple bisection search
            if (scalar < controlPoints[guess].scalar)
            {
                high = guess;
            } else {
                low = guess;
            }
        }
        return low;
    }
    
    void copyColor(float to[3], const float from[3])
    {
        for (int i = 0; i < 3; ++i)
        {
            CaretAssert(from[i] >= 0.0f && from[i] <= 1.0f);
            to[i] = from[i];
        }
    }
}

PaletteNew::PaletteNew(vector<ScalarColor> posRange, float zeroColor[3], vector<ScalarColor> negRange) : m_posRange(posRange), m_negRange(negRange)
{
    CaretAssert(posRange[0].scalar == 0.0f);
    CaretAssert(posRange.back().scalar == 1.0f);
    CaretAssert(negRange[0].scalar == -1.0f);
    CaretAssert(negRange.back().scalar == 0.0f);
    copyColor(m_zeroColor, zeroColor);
}

void PaletteNew::getPositiveColor(const float scalar, float rgbOut[3]) const
{
    CaretAssert(scalar >= 0.0f);
    m_posRange.getPaletteColor(scalar, rgbOut);
}

void PaletteNew::getNegativeColor(const float scalar, float rgbOut[3]) const
{
    CaretAssert(scalar <= 0.0f);
    m_negRange.getPaletteColor(scalar, rgbOut);
}

PaletteNew::PaletteRange::PaletteRange(const vector<ScalarColor> controlPoints)
{
    CaretAssert(controlPoints.size() > 1);
    for (size_t i = 1; i < controlPoints.size(); ++i)
    {
        if (controlPoints[i].scalar <= controlPoints[i - 1].scalar)
        {
            CaretAssert(false);
            throw CaretException("palette ranges must use unique, ascending-sorted scalars");
        }
    }
    m_controlPoints = controlPoints;
    m_lowPoint = m_controlPoints[0].scalar;
    m_highPoint = m_controlPoints.back().scalar;
    float diff = m_highPoint - m_lowPoint;//in practice this should always be 1
    m_lookup[0] = 0;
    m_lookup[BUCKETS] = m_controlPoints.size() - 2;//do not allow search to ever say that it is past the last control point
    for (int i = 1; i < BUCKETS; ++i)
    {
        m_lookup[i] = searchLowSide(m_lowPoint + (diff * i) / (BUCKETS), m_controlPoints, 0, m_controlPoints.size() - 1);//save the low control point at each bucket wall (can derive high from next bucket)
    }
}

void PaletteNew::PaletteRange::getPaletteColor(const float scalar, float rgbOut[3]) const
{
    CaretAssert(MathFunctions::isNumeric(scalar));
    int bucket = int((scalar - m_lowPoint) / (m_highPoint - m_lowPoint) * BUCKETS);//this implicitly tests for below or above the palette range
    if (bucket < 0)
    {
        copyColor(rgbOut, m_controlPoints[0].color);
        return;
    }
    if (bucket >= BUCKETS)
    {
        copyColor(rgbOut, m_controlPoints.back().color);
        return;
    }
    int lowSide = searchLowSide(scalar, m_controlPoints, m_lookup[bucket], m_lookup[bucket + 1] + 1);//this is why m_lookup is [BUCKETS + 1]
    float interpFactor = (scalar - m_controlPoints[lowSide].scalar) / (m_controlPoints[lowSide + 1].scalar - m_controlPoints[lowSide].scalar);
    for (int i = 0; i < 3; ++i)
    {//never output outside of [0, 1]
        rgbOut[i] = max(0.0f, min(1.0f, m_controlPoints[lowSide].color[i] * (1 - interpFactor) + m_controlPoints[lowSide + 1].color[i] * interpFactor));
    }
}

vector<float> PaletteNew::PaletteRange::getPerceptualGradient(const int numBuckets) const
{
    CaretAssert(numBuckets >= 3);//one at min, one at max, all equally spaced
    vector<float> ret(numBuckets);
    ScalarColor lowColor, midColor, highColor;
    midColor.scalar = m_lowPoint;
    getPaletteColor(midColor.scalar, midColor.color);
    lowColor = midColor;//use central difference, with unsymmetric endpoints
    for (int i = 1; i < numBuckets - 1; ++i)
    {
        float interpFactor = float(i) / (numBuckets - 1);
        highColor.scalar = (1 - interpFactor) * m_lowPoint + interpFactor * m_highPoint;
        getPaletteColor(highColor.scalar, highColor.color);
        ret[i - 1] = ColorFunctions::perceptualDistanceSRGB(highColor.color, lowColor.color) / (highColor.scalar - lowColor.scalar);
        lowColor = midColor;
        midColor = highColor;
    }
    ret[numBuckets - 1] = ColorFunctions::perceptualDistanceSRGB(highColor.color, midColor.color) / (highColor.scalar - midColor.scalar);
    return ret;
}
