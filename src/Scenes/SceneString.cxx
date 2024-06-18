
/*LICENSE_START*/
/*
 *  Copyright (C) 2014  Washington University School of Medicine
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

#include <limits>

#define __SCENE_STRING_DECLARE__
#include "SceneString.h"
#undef __SCENE_STRING_DECLARE__

using namespace caret;


    
/**
 * \class caret::SceneString 
 * \brief For storage of a string value in a scene.
 * \ingroup Scene
 *
 * See the documentation in the class Scene for how to use the Scene system.
 */

/**
 * Constructor.
 *
 * @param name
 *   Name of object.
 * @param value
 *   Value of object.
 */
SceneString::SceneString(const AString& name,
                           const AString& value)
: ScenePrimitive(name,
                 SceneObjectDataTypeEnum::SCENE_STRING)
{
    m_value = value;
}

SceneString::SceneString(const SceneString& rhs): ScenePrimitive(rhs.getName(), SceneObjectDataTypeEnum::SCENE_STRING)
{
    m_value = rhs.m_value;
}

SceneObject* SceneString::clone() const
{
    return new SceneString(*this);
}

/**
 * Destructor.
 */
SceneString::~SceneString()
{
    
}

/**
 * Set the value.
 * @param value
 *   New value.
 */
void 
SceneString::setValue(const AString& value)
{
    m_value = value;
}

/**
 * @return The value as a boolean data type.
 */
bool 
SceneString::booleanValue() const
{
    m_restoredFlag = true;
    const bool b = m_value.toBool();
    return b;
}

/**
 * @return The value as a float data type.
 * If the string does not convert to a float number,
 * 0.0 is returned.
 */
float
SceneString::floatValue() const
{
    m_restoredFlag = true;
    bool isValid = false;
    float f = m_value.toFloat(&isValid);
    if (isValid == false) {
        f = 0.0;
    }
    return f;
}

/**
 * @return The value as a integer data type.
 * If the string does not convert to an integer number,
 * 0 is returned.
 */
int32_t 
SceneString::integerValue() const
{
    m_restoredFlag = true;
    bool isValid = false;
    int32_t i = m_value.toInt(&isValid);
    if (isValid == false) {
        i = 0;
    }
    return i;
}

/**
 * @return The value as a long integer data type.
 */
int64_t
SceneString::longIntegerValue() const
{
    m_restoredFlag = true;
    bool isValid = false;
    int64_t i = m_value.toLong(&isValid);
    if (isValid == false) {
        i = 0;
    }
    return i;
}

/**
 * @return The value as a string data type.
 */
AString 
SceneString::stringValue() const
{
    m_restoredFlag = true;
    return m_value;
}

/**
 * Get the values as an unsigned byte.
 * @param arrayIndex
 *    Index of element.
 * @return The value.
 */
uint8_t
SceneString::unsignedByteValue() const
{
    m_restoredFlag = true;
    bool isValid = false;
    const uint32_t i = m_value.toUInt(&isValid);
    if ( ! isValid) {
        return 0;
    }
    
    /*
     * Note since "i" is unsigned there is no need
     * to compare with std::numeric_limits<uint8_t>::min()
     * since it is zero.
     */
    if (i > std::numeric_limits<uint8_t>::max()) {
        return std::numeric_limits<uint8_t>::max();
    }
    
    const uint8_t b = static_cast<uint8_t>(i);
    return b;
}

