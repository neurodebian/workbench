
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

#define __SCENE_LONG_INTEGER_DECLARE__
#include "SceneLongInteger.h"
#undef __SCENE_LONG_INTEGER_DECLARE__

using namespace caret;


    
/**
 * \class caret::SceneLongInteger
 * \brief For storage of a longinteger value in a scene.
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
SceneLongInteger::SceneLongInteger(const AString& name,
                       const int64_t value)
: ScenePrimitive(name,
                 SceneObjectDataTypeEnum::SCENE_INTEGER)
{
    m_value = value;
}

SceneLongInteger::SceneLongInteger(const SceneLongInteger& rhs): ScenePrimitive(rhs.getName(), SceneObjectDataTypeEnum::SCENE_INTEGER)
{
    m_value = rhs.m_value;
}

SceneObject* SceneLongInteger::clone() const
{
    return new SceneLongInteger(*this);
}

/**
 * Destructor.
 */
SceneLongInteger::~SceneLongInteger()
{
    
}

/**
 * Set the value.
 * @param value
 *    The new value.
 */
void 
SceneLongInteger::setValue(const int64_t value)
{
    m_value = value;
}

/**
 * @return The value as a boolean data type.
 */
bool 
SceneLongInteger::booleanValue() const
{
    m_restoredFlag = true;
    const bool b = ((m_value != 0.0) ? true : false);
    return b;
}

/**
 * @return The value as a float data type.
 */
float
SceneLongInteger::floatValue() const
{
    m_restoredFlag = true;
    const float f = static_cast<float>(m_value);
    return f;
}

/**
 * @return The value as a integer data type.
 */
int32_t 
SceneLongInteger::integerValue() const
{
    m_restoredFlag = true;
    return m_value;
}

/**
 * @return The value as a long integer data type.
 */
int64_t
SceneLongInteger::longIntegerValue() const
{
    m_restoredFlag = true;
    return m_value;
}

/**
 * @return The value as a string data type.
 */
AString 
SceneLongInteger::stringValue() const
{
    m_restoredFlag = true;
    const AString s = AString::number(m_value);
    return s;
}

/**
 * Get the values as an unsigned byte.
 * @param arrayIndex
 *    Index of element.
 * @return The value.
 */
uint8_t
SceneLongInteger::unsignedByteValue() const
{
    m_restoredFlag = true;
    if (m_value > std::numeric_limits<uint8_t>::max()) {
        return std::numeric_limits<uint8_t>::max();
    }
    else if (m_value < std::numeric_limits<uint8_t>::min()) {
        return std::numeric_limits<uint8_t>::min();
    }
    
    const uint8_t b = static_cast<uint8_t>(m_value);
    return b;
}
