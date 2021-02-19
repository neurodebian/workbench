
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

#define __FILE_IDENTIFICATION_ATTRIBUTES_DECLARE__
#include "FileIdentificationAttributes.h"
#undef __FILE_IDENTIFICATION_ATTRIBUTES_DECLARE__

#include "CaretAssert.h"
#include "SceneClass.h"
#include "SceneClassAssistant.h"

using namespace caret;


    
/**
 * \class caret::FileIdentificationAttributes 
 * \brief File identification attributes used in identification file filtering
 * \ingroup Files
 */

/**
 * Constructor.
 */
FileIdentificationAttributes::FileIdentificationAttributes()
: CaretObject()
{
    reset();
    
    m_sceneAssistant = std::unique_ptr<SceneClassAssistant>(new SceneClassAssistant());
    m_sceneAssistant->add("m_enabled",
                          &m_enabled);
    m_sceneAssistant->add<FileIdentificationMapSelectionEnum, FileIdentificationMapSelectionEnum::Enum>("m_mapSelectionMode",
                                                                                                        &m_mapSelectionMode);
    m_sceneAssistant->add("m_mapIndex",
                          &m_mapIndex);
}

void
FileIdentificationAttributes::reset()
{
    m_enabled = false;
    m_mapSelectionMode = FileIdentificationMapSelectionEnum::SELECTED;
    m_mapIndex = 0;
}

/**
 * Destructor.
 */
FileIdentificationAttributes::~FileIdentificationAttributes()
{
}

/**
 * Copy constructor.
 * @param obj
 *    Object that is copied.
 */
FileIdentificationAttributes::FileIdentificationAttributes(const FileIdentificationAttributes& obj)
: CaretObject(obj),
SceneableInterface(obj)
{
    this->copyHelperFileIdentificationAttributes(obj);
}

/**
 * Assignment operator.
 * @param obj
 *    Data copied from obj to this.
 * @return 
 *    Reference to this object.
 */
FileIdentificationAttributes&
FileIdentificationAttributes::operator=(const FileIdentificationAttributes& obj)
{
    if (this != &obj) {
        CaretObject::operator=(obj);
        this->copyHelperFileIdentificationAttributes(obj);
    }
    return *this;    
}

/**
 * Helps with copying an object of this type.
 * @param obj
 *    Object that is copied.
 */
void 
FileIdentificationAttributes::copyHelperFileIdentificationAttributes(const FileIdentificationAttributes& obj)
{
    m_enabled = obj.m_enabled;
    m_mapSelectionMode = obj.m_mapSelectionMode;
    m_mapIndex = obj.m_mapIndex;
}

/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
FileIdentificationAttributes::toString() const
{
    return "FileIdentificationAttributes";
}

/**
 * @return enabled for identification
 */
bool
FileIdentificationAttributes::isEnabled() const
{
    return m_enabled;
}

/**
 * Set enabled for identification
 *
 * @param enabled
 *    New value for enabled for identification
 */
void
FileIdentificationAttributes::setEnabled(const bool enabled)
{
    m_enabled = enabled;
}

/**
 * @return Map selection mode
 */
FileIdentificationMapSelectionEnum::Enum
FileIdentificationAttributes::getMapSelectionMode() const
{
    return m_mapSelectionMode;
}

/**
 * Set the map selection mode
 *
 * @param mapSelectionMode
 *    New map selection mode
 */
void
FileIdentificationAttributes::setMapSelectionMode(const FileIdentificationMapSelectionEnum::Enum mapSelectionMode)
{
    m_mapSelectionMode = mapSelectionMode;
}

/**
 * @return map selected for identification
 */
int32_t
FileIdentificationAttributes::getMapIndex() const
{
    return m_mapIndex;
}

/**
 * Set map selected for identification
 *
 * @param mapIndex
 *    New value for map selected for identification
 */
void
FileIdentificationAttributes::setMapIndex(const int32_t mapIndex)
{
    m_mapIndex = mapIndex;
}

/**
 * Save information specific to this type of model to the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    saving the scene.
 *
 * @param instanceName
 *    Name of instance in the scene.
 */
SceneClass*
FileIdentificationAttributes::saveToScene(const SceneAttributes* sceneAttributes,
                                 const AString& instanceName)
{
    SceneClass* sceneClass = new SceneClass(instanceName,
                                            "FileIdentificationAttributes",
                                            1);
    m_sceneAssistant->saveMembers(sceneAttributes,
                                  sceneClass);
    
    // Uncomment if sub-classes must save to scene
    //saveSubClassDataToScene(sceneAttributes,
    //                        sceneClass);
    
    return sceneClass;
}

/**
 * Restore information specific to the type of model from the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass from which model specific information is obtained.
 */
void
FileIdentificationAttributes::restoreFromScene(const SceneAttributes* sceneAttributes,
                                      const SceneClass* sceneClass)
{
    reset();
    
    if (sceneClass == NULL) {
        return;
    }
    
    m_sceneAssistant->restoreMembers(sceneAttributes,
                                     sceneClass);    
    
    //Uncomment if sub-classes must restore from scene
    //restoreSubClassDataFromScene(sceneAttributes,
    //                             sceneClass);
    
}

