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

/*
 * When GLEW is used, CaretOpenGLInclude.h will include "Gl/glew.h".
 * Gl/glew.h MUST BE BEFORE Gl/gl.h and Gl/gl.h is included by
 * some Qt classes so, we must include CaretOpenGL.h first.
 */
#include "CaretOpenGLInclude.h"

#include <QApplication>
#include <QDesktopWidget>
#include <QLabel>
#ifndef WORKBENCH_USE_QT5_QOPENGL_WIDGET
#include <QGLPixelBuffer>
#endif
#include <QSplashScreen>
#include <QStyleFactory>
#include <QThread>

#include <cstdlib>
#include <ctime>
#include <iostream>

#include "ApplicationInformation.h"
#include "BrainBrowserWindow.h"
#include "BrainOpenGL.h"
#include "BrainOpenGLWidget.h"
#include "CaretAssert.h"
#include "CaretCommandLine.h"
#include "CaretHttpManager.h"
#include "CaretLogger.h"
#include "CaretPreferences.h"
#include "CommandOperationManager.h"
#include "EventBrowserWindowNew.h"
#include "EventManager.h"
#include "FileInformation.h"
#include "GuiManager.h"
#include "MacApplication.h"
#include "ProgramParameters.h"
#include "RecentFilesDialog.h"
#include "RecentFilesSystemAccessModeEnum.h"
#include "SceneDialog.h"
#include "SessionManager.h"
#include "SystemUtilities.h"
#include "WorkbenchQtMessageHandler.h"
#include "WuQMessageBox.h"
#include "WuQtUtilities.h"

using namespace caret;
using namespace std;


//struct for communicating stuff back to main from parseCommandLine
struct ProgramState
{
    vector<AString> fileList;
    int windowSizeXY[2];
    int windowPosXY[2];
    int graphicsSizeXY[2];
    bool showSplash;
    bool macMenuFlag = false;
    
    AString directoryName;
    
    AString specFileNameLoadWithDialog;
    AString specFileNameLoadAll;
    
    AString sceneFileName;
    AString sceneFileNameNoDialog;
    AString sceneNameOrNumber;
    
    ProgramState();
};

class GuiBeeper : public WorkbenchQtMessageHandler::Beeper
{
public:
    GuiBeeper() { }
    
    ~GuiBeeper() { }
    
    /**
     * Override to make the beep sound
     */
    virtual void makeBeep() override { GuiManager::get()->beep(); }
    
};

//declare the functions associated with command line
void printHelp(const AString& progName);
void parseCommandLine(const AString& progName, ProgramParameters* myParams, ProgramState& myState);
int 
main(int argc, char* argv[])
{
    srand(time(NULL));
    int result;
    {
        /*
        * Handle uncaught exceptions
        */
        SystemUtilities::setUnexpectedHandler();
        
        /*
        * Create the session manager.
        */
        SessionManager::createSessionManager(ApplicationTypeEnum::APPLICATION_TYPE_GRAPHICAL_USER_INTERFACE);

        /*
        * Parameters for the program.
        */
        ProgramParameters parameters(argc, argv);
        caret_global_commandLine_init(parameters);

        //begin parsing command line
        ProgramState myState;
        FileInformation progInfo(argv[0]);
        AString progName = progInfo.getFileName();
        parseCommandLine(progName, &parameters, myState);
        
        /*
        * Log the command parameters.
        */
        CaretLogFine("Running: " + caret_global_commandLine);
        
        /*
         * Setup OpenGL if using Qt 5's QOpenGLWidget.
         * QOpenGLWidget's documentation indicates that
         * default format must be initizlied before
         * QApplication is created.
         */
#ifdef WORKBENCH_USE_QT5_QOPENGL_WIDGET
        QCoreApplication::setAttribute(Qt::AA_ShareOpenGLContexts);
        BrainOpenGLWidget::initializeDefaultGLFormat();
#endif
        
        //change the default graphics system on mac to avoid rendering performance issues with qwtplotter
#ifdef CARET_OS_MACOSX
        // Qt-Deprecated: QApplication::setGraphicsSystem("raster");
        MacApplication app(argc, argv);
#else //CARET_OS_MACOSX        
        QApplication app(argc, argv);
#endif //CARET_OS_MACOSX

        /*
         * Create the GUI Manager.
         * Moved here as part of WB-842.  In OSX Mojave (10.14),
         * the open file event appears to be delivered very quickly
         * so the GUI manager need to be created sooner than
         * before.
         */
        GuiManager::createGuiManager();
        app.processEvents();

        ApplicationInformation applicationInformation;

        QApplication::addLibraryPath(
            QApplication::applicationDirPath()
            + QDir::separator()
            + "plugins");
        QApplication::setApplicationName(applicationInformation.getName());
        QApplication::setApplicationVersion(applicationInformation.getVersion());
        QApplication::setOrganizationDomain("brainvis.wustl.edu");
        QApplication::setOrganizationName("Van Essen Lab");
        /*
         * Override the system local to US - English
         */
        QLocale::setDefault(QLocale(QLocale::English, QLocale::UnitedStates));
        
        /*
        * Make sure OpenGL is available.
        */
        const QString noOpenGLMessage("OpenGL (3D Graphics System) is not available.  "
                                      "This may be caused by missing or outdated OpenGL libraries."
        
#ifdef CARET_OS_LINUX
                                      "On Linux, this may be caused by a missing plugin "
                                      "library <plugins-path>/xcbglintegrations/libqxcb-glx-integration.so."
#endif
                                      );
        
                                      
#ifndef WORKBENCH_USE_QT5_QOPENGL_WIDGET
        if ( ! QGLFormat::hasOpenGL()) {
            app.processEvents();
            WuQMessageBox::errorOk(NULL,
                                   noOpenGLMessage);
            app.processEvents();
            
            return -1;
        }
#endif
            
        /*
        * Setup OpenGL if NOT using Qt 5's QOpenGLWidget
        */
#ifndef WORKBENCH_USE_QT5_QOPENGL_WIDGET
        BrainOpenGLWidget::initializeDefaultGLFormat();
#endif
        
        /*
         * Setup handlers (callbacks) for messages produced by Qt
         */
        GuiBeeper beeper;
        WorkbenchQtMessageHandler::setupHandler(&beeper);

        /*
         * Log debug status
         */
        CaretLogConfig("Version " + applicationInformation.getVersion());
        CaretLogConfig(applicationInformation.getCompiledWithDebugStatus());


        /*
         * Mac Option to put menus on window and not use the native tool bar
         */
        if (myState.macMenuFlag) {
            app.setAttribute(Qt::AA_DontUseNativeMenuBar);
        }
        
        /*
         * Enabled the desired splash screen based upon user preferences
         * and command line options.  Do not show selection splash screen
         * if there has listed files on the command line.
         */ 
        CaretPreferences* preferences = SessionManager::get()->getCaretPreferences();
        bool showSelectionSplashScreen = preferences->isSplashScreenEnabled();
        bool showImageSplashScreen = (! showSelectionSplashScreen);
        if ( ! myState.showSplash) {
            showSelectionSplashScreen = false;
            showImageSplashScreen = false;
        }
        
        /*
         * DISABLE IMAGE SPLASH SCREEN
         */
        showImageSplashScreen = false;
        
        /*
        * Splash Screen
        */
        QPixmap splashPixmap;
        QSplashScreen splashScreen;
        if (showImageSplashScreen) {
            if (WuQtUtilities::loadPixmap(":/Splash/hcp.png", splashPixmap)) {
                splashScreen.setPixmap(splashPixmap);
                splashScreen.showMessage("Starting Workbench...");
                splashScreen.show();
                app.processEvents();
                SystemUtilities::sleepSeconds(2);
            }
        }
        
        /*
         * Letting the App process events will allow the message for a 
         * double-clicked spec file in Mac OSX to get processed.
         */
        app.processEvents();
        
        /*
         * Now that events have processed, see if there was a request for
         * a data file to open.
         */
        const AString dataFileNameFromOS = GuiManager::get()->getNameOfDataFileToOpenAfterStartup();
        if ( ! dataFileNameFromOS.isEmpty()) {
            showSelectionSplashScreen = false;
            if (dataFileNameFromOS.endsWith(DataFileTypeEnum::toFileExtension(DataFileTypeEnum::SPECIFICATION))) {
                myState.specFileNameLoadWithDialog = dataFileNameFromOS;
            }
            else {
                myState.fileList.push_back(dataFileNameFromOS);
            }
        }
        
        /*
         * Show file selection splash screen if enabled via user's preferences
         */
        if (showSelectionSplashScreen) {
            AString dataFileName;
            int32_t sceneIndex(-1);
            RecentFilesDialog::ResultModeEnum result = RecentFilesDialog::runDialog(RecentFilesDialog::RunMode::SPLASH_SCREEN,
                                                                                    dataFileName,
                                                                                    sceneIndex);
            switch (result) {
                case RecentFilesDialog::ResultModeEnum::CANCEL:
                    break;
                case RecentFilesDialog::ResultModeEnum::LOAD_FILES_IN_SPEC_FILE:
                    myState.specFileNameLoadAll = dataFileName;
                    break;
                case RecentFilesDialog::ResultModeEnum::LOAD_SCENE_FROM_SCENE_FILE:
                    myState.sceneFileName     = dataFileName;
                    myState.sceneNameOrNumber = AString::number(sceneIndex);
                    break;
                case RecentFilesDialog::ResultModeEnum::OPEN_DIRECTORY:
                    myState.directoryName = dataFileName;
                    break;
                case RecentFilesDialog::ResultModeEnum::OPEN_FILE:
                    if ( ! dataFileName.isEmpty()) {
                        bool validFlag(false);
                        if (DataFileTypeEnum::fromFileExtension(dataFileName,
                                                                &validFlag) == DataFileTypeEnum::SPECIFICATION) {
                            myState.specFileNameLoadWithDialog = dataFileName;
                        }
                        else {
                            myState.fileList.push_back(dataFileName);
                        }
                    }
                    break;
                case RecentFilesDialog::ResultModeEnum::OPEN_OTHER:
                    myState.directoryName = dataFileName;
                    break;
            }
        }
        
        /*
        * Create and display a main window.
        * If not done as pointer, the event object is listed as an
        * object that was not deleted by CaretObject::printListOfObjectsNotDeleted
        * since it does not go out of scope.
        */
        
        EventBrowserWindowNew newBrowserWindow(NULL, NULL);
        EventManager::get()->sendEvent(newBrowserWindow.getPointer());
        splashScreen.close();

        
        BrainBrowserWindow* myWindow = GuiManager::get()->getBrowserWindowByWindowIndex(0);
        if ( ! myWindow->hasValidOpenGL()) {
            app.processEvents();
            WuQMessageBox::errorOk(NULL,
                                   noOpenGLMessage);
            app.processEvents();
            
            return -1;
        }
        
        /**
         * Test for the required OpenGL version is available.
         */
        AString requiredOpenGLMessage;
        if ( ! BrainOpenGL::testForRequiredOpenGLVersion(requiredOpenGLMessage)) {
            app.processEvents();
            if ( ! WuQMessageBox::warningAcceptReject(NULL,
                                                      requiredOpenGLMessage,
                                                      "Continue",
                                                      "Exit")) {
                return -1;
            }
        }
        
        /*
         * Initialize GLEW if it is being used
         */
        AString glewErrorMessage;
#ifdef HAVE_GLEW
        GLenum err = glewInit();
        if (GLEW_OK != err) {
            /* Problem: glewInit failed, something is seriously wrong. */
            
            glewErrorMessage = ("GLEW failed to initialize.\n"
                                + AString(reinterpret_cast<const char*>(glewGetErrorString(err))));
            CaretLogSevere(glewErrorMessage);
        }
        CaretLogInfo("Status: Using GLEW version "
                     + AString(reinterpret_cast<const char*>(glewGetString(GLEW_VERSION))));
#endif /* HAVE_GLEW */
        if ( ! glewErrorMessage.isEmpty()) {
            app.processEvents();
            WuQMessageBox::errorOk(NULL,
                                   glewErrorMessage);
            app.processEvents();
            
            return -1;
        }
        
        if ((myState.windowSizeXY[0] > 0) && (myState.windowSizeXY[1] > 0))
        {
            if (myState.windowPosXY[0] > 0 && myState.windowPosXY[1] > 0)
            {
                myWindow->setGeometry(myState.windowPosXY[0], myState.windowPosXY[1],
                                      myState.windowSizeXY[0], myState.windowSizeXY[1]);
            } else {
                myWindow->setFixedSize(myState.windowSizeXY[0],
                                       myState.windowSizeXY[1]);
            }
        } else {
            if (myState.windowPosXY[0] > 0 && myState.windowPosXY[1] > 0)
            {
                myState.windowSizeXY[0] = myWindow->width();
                myState.windowSizeXY[1] = myWindow->height();
                myWindow->setGeometry(myState.windowPosXY[0], myState.windowPosXY[1],
                                      myState.windowSizeXY[0], myState.windowSizeXY[1]);
            }
        }
        
        if (myState.graphicsSizeXY[0] > 0 && myState.graphicsSizeXY[1] > 0) {
            myWindow->setGraphicsWidgetFixedSize(myState.graphicsSizeXY[0], myState.graphicsSizeXY[1]);
        }
        
        if ( ! myState.specFileNameLoadAll.isEmpty()) {
            myWindow->loadFilesFromCommandLine({ myState.specFileNameLoadAll },
                                               BrainBrowserWindow::LOAD_SPEC_FILE_CONTENTS_VIA_COMMAND_LINE);
        }
        else if ( ! myState.specFileNameLoadWithDialog.isEmpty()) {
            myWindow->loadFilesFromCommandLine({ myState.specFileNameLoadWithDialog },
                                               BrainBrowserWindow::LOAD_SPEC_FILE_WITH_DIALOG_VIA_COMMAND_LINE);
        }
        
        if (! myState.fileList.empty()) {
            myWindow->loadFilesFromCommandLine(myState.fileList,
                                               BrainBrowserWindow::LOAD_SPEC_FILE_WITH_DIALOG);
        }
        
        if ( ! myState.sceneFileName.isEmpty()) {
            myWindow->loadSceneFromCommandLine(myState.sceneFileName,
                                               myState.sceneNameOrNumber,
                                               BrainBrowserWindow::LoadSceneFromCommandLineDialogMode::SHOW_YES);
        }
        
        if ( ! myState.sceneFileNameNoDialog.isEmpty()) {
            myWindow->loadSceneFromCommandLine(myState.sceneFileNameNoDialog,
                                               myState.sceneNameOrNumber,
                                               BrainBrowserWindow::LoadSceneFromCommandLineDialogMode::SHOW_NO);
        }
        
        if ( ! myState.directoryName.isEmpty()) {
            myWindow->loadDirectoryFromCommandLine(myState.directoryName);
        }
        
#ifndef WORKBENCH_USE_QT5_QOPENGL_WIDGET
        if (QGLPixelBuffer::hasOpenGLPbuffers()) {
            CaretLogConfig("OpenGL PBuffers are supported");
        }
        else {
            CaretLogConfig("OpenGL PBuffers are NOT supported");
        }
#endif
    
        /*
         * Log local (language, country)
         */
        QLocale sytemLocale = QLocale::system();
        CaretLogConfig("Local Language="
                       + QLocale::languageToString(sytemLocale.language())
                       + " Country="
                       + QLocale::countryToString(sytemLocale.country()));

        /*
         * Resolution of screens
         */
        AString screenSizeText = "Screen Sizes: ";
        QDesktopWidget* dw = QApplication::desktop();
        const int numScreens = dw->screenCount();
        for (int i = 0; i < numScreens; i++) {
            const QRect rect = dw->screenGeometry(i);
            const int x = rect.x();
            const int y = rect.y();
            const int w = rect.width();
            const int h = rect.height();
            screenSizeText.appendWithNewLine("Screen index="
                                             + AString::number(i)
                                             + ", x="
                                             + AString::number(x)
                                             + ", y="
                                             + AString::number(y)
                                             + ", w="
                                             + AString::number(w)
                                             + ", h="
                                             + AString::number(h));
        }
        screenSizeText.appendWithNewLine("Primary Screen="
                                         + AString::number(dw->primaryScreen()));
        if (dw->isVirtualDesktop()) {
            screenSizeText.appendWithNewLine("Virtual Desktop=YES");
        }
        else {
            screenSizeText.appendWithNewLine("Virtual Desktop=NO");
        }
        QWidget* screenWidget = dw->screen();
        QRect screenWidgetRect = screenWidget->geometry();
        screenSizeText.appendWithNewLine("Desktop: x="
                                         + AString::number(screenWidgetRect.x())
                                         + ", y="
                                         + AString::number(screenWidgetRect.y())
                                         + ", w="
                                         + AString::number(screenWidgetRect.width())
                                         + ", h="
                                         + AString::number(screenWidgetRect.height()));
        
        screenSizeText.appendWithNewLine("Logical DPI: x="
                                         + AString::number(dw->logicalDpiX())
                                         + ", y="
                                         + AString::number(dw->logicalDpiY()));
        
        screenSizeText.appendWithNewLine("Physical DPI: x="
                                         + AString::number(dw->physicalDpiX())
                                         + ", y="
                                         + AString::number(dw->physicalDpiY()));
        
        screenSizeText.appendWithNewLine("Width/height (mm): x="
                                         + AString::number(dw->widthMM())
                                         + ", y="
                                         + AString::number(dw->heightMM()));
        
        CaretLogConfig(screenSizeText);
        
        /*
        * Start the app which will launch the main window.
        */
        result = app.exec();
        
        /*
        * Hiding the window removes it from the event loop on Windows, which is necessary to
        * prevent paint events from causing assertion errors when the Window is destroyed
        * Although this is a Window's only bug, it's probably good practice to do on all platforms
        */
        //theMainWindow->hide();
        
        /*
        * Delete the GUI Manager.
        */
        GuiManager::deleteGuiManager();
        
        /*
         * Delete the command manager
         */
        CommandOperationManager::deleteCommandOperationManager();
        
        /*
        * Delete the session manager.
        */
        SessionManager::deleteSessionManager();
        
        CaretHttpManager::deleteHttpManager();
    }
    
    /*
     * See if any objects were not deleted.
     */
    CaretObject::printListOfObjectsNotDeleted(true);
    //FreeConsole();
    
    return result;
}

void printHelp(const AString& progName)
{
    cout
    << "Usage: " << progName << " [options] [files]" << endl
    << endl
    << "    [files], if present, can be a single spec file, or multiple data files" << endl
    << endl
    << "Options:" << endl
    << "    -help" << endl
    << "        display this usage text" << endl
    << endl
    << "    -graphics-size  <X Y>" << endl
    << "        Set the size of the graphics region." << endl
    << "        If this option is used you WILL NOT be able" << endl
    << "        to change the size of the graphic region. It" << endl
    << "        may be useful when image captures of a particular" << endl
    << "        size are desired." << endl
    << endl
    << "    -logging <level>" << endl
    << "       Set the logging level." << endl
    << "       Valid Levels are:" << endl;
    
    std::vector<LogLevelEnum::Enum> logLevels;
    LogLevelEnum::getAllEnums(logLevels);
    for (std::vector<LogLevelEnum::Enum>::iterator iter = logLevels.begin();
         iter != logLevels.end();
         iter++) {
        cout << "           " << qPrintable(LogLevelEnum::toName(*iter)) << endl;
    }
    
    cout
    << endl
    << "    -mac-menu-in-window" << endl
    << "        MacOS Only - Menus are at the top of each main window  "
    << "        similar to Linux and Window Applications.  "
    << "        May be useful for creating tutorial images and videos." << endl
    << endl
    << "    -no-recent-files-dialog" << endl
    << "         Inhibits display of Open Recent Files Dialog at wb_view startup." << endl
    << endl
    << "    -no-splash" << endl
    << "        (Obsolete) Replaced by \"-no-recent-files-dialog\"" << endl
    << "        Splash screen was replaced with Open Recent Fiels Dialog." << endl
    << endl
    << "    -recent-files-mode <mode>" << endl
    << "        Set the recent file's file system access mode" << endl
    << "        (overrides and replaces value in preferences)." << endl
    << "        Using an off mode prevents file system access for " << endl
    << "        obtaining last modified time and file existance." << endl
    << "        This option may be useful if recent files are on a mounted" << endl
    << "        file system that is having problems that may cause wb_view" << endl
    << "        to hang at startup.  Valid modes are:" << endl;
    
    std::vector<RecentFilesSystemAccessModeEnum::Enum> recentFilesModes;
    RecentFilesSystemAccessModeEnum::getAllEnums(recentFilesModes);
    for (auto rfm : recentFilesModes) {
        std::cout  << "            " <<RecentFilesSystemAccessModeEnum::toName(rfm) << endl;
    }
    
    cout
    << endl
    << "    -scene-load <scene-file-name> <scene-name-or-number>" << endl
    << "        load the specified scene file and display the scene " << endl
    << "        in the file that matches by name or number.  Name" << endl
    << "        takes precedence over number.  The scene numbers " << endl
    << "        start at one.  The scene dialog remains visible " << endl
    << "        after loading of the scene." << endl
    << endl
    << "    -scene-load-hd <scene-file-name> <scene-name-of-number>" << endl
    << "        Same as \"-scene-load\" except that the scene dialog " << endl
    << "        is hidden after the scene has loaded." << endl
    << endl
    << "    -style <style-name>" << endl
    << "        change the window style to the specified style" << endl
    << "        the following styles are valid on this system:" << endl;
    
    QStringList styleList = QStyleFactory::keys();
    QStringListIterator styleListIterator(styleList);
    while (styleListIterator.hasNext()) {
        cout << "           " << qPrintable(styleListIterator.next()) << endl;
    }
    
    cout
    << "        The selected style is listed on the About wb_view dialog" << endl
    << "        available from the File Menu (On Macs: wb_view Menu). " << endl
    << "        Press the \"More\" button to see the selected style." << endl
    << "        Other styles may be available on other systems." << endl
    << endl
    << "    -spec-load-all" << endl
    << "        load all files in the given spec file, don't show spec file dialog" << endl
    << endl
    << "    -window-size  <X Y>" << endl
    << "        Set the size of the browser window" << endl
    << endl
    << "    -window-pos  <X Y>" << endl
    << "        Set the position of the browser window" << endl
    << endl;
}

void parseCommandLine(const AString& progName, ProgramParameters* myParams, ProgramState& myState)
{
    bool hasFatalError = false;
    
    const AString moreThanOneSpecFileErrorMessage("More than one spec file is NOT allowed in options");
    
    try {
        while (myParams->hasNext())
        {
            AString thisParam = myParams->nextString("option");
            if (thisParam[0] == '-')
            {
                if (thisParam == "-style")
                {
                    myParams->nextString("style");//discard, QApplication handles this
                } else if (thisParam == "-help") {
                    printHelp(progName);
                    exit(0);
                } else if (thisParam == "-logging") {
                    if (myParams->hasNext()) {
                        const AString logLevelName = myParams->nextString("Logging Level").toUpper();
                        bool valid = false;
                        const LogLevelEnum::Enum level = LogLevelEnum::fromName(logLevelName, &valid);
                        if (valid)
                        {
                            /*
                             * Note settings logging level in preferences will also
                             * set logging level in the caret logger.
                             */
                            CaretPreferences* prefs = SessionManager::get()->getCaretPreferences();
                            prefs->setLoggingLevel(level);
                        }
                        else {
                            cerr << "Invalid logging level \""
                            << qPrintable(logLevelName)
                            << "\" for \"-logging\" option" << std::endl;
                            hasFatalError = true;
                        }
                    }
                } else if (thisParam == "-mac-menu-in-window") {
                    myState.macMenuFlag = true;
                } else if ((thisParam == "-no-splash")
                           || (thisParam == "-no-recent-files-dialog")) {
                    myState.showSplash = false;
                } else if (thisParam == "-recent-files-mode") {
                    if (myParams->hasNext()) {
                        const AString recentFilesModeName = myParams->nextString("Recent Files Mode").toUpper();
                        bool valid = false;
                        const RecentFilesSystemAccessModeEnum::Enum recentFilesMode = RecentFilesSystemAccessModeEnum::fromName(recentFilesModeName, &valid);
                        if (valid)
                        {
                            CaretPreferences* prefs = SessionManager::get()->getCaretPreferences();
                            prefs->setRecentFilesSystemAccessMode(recentFilesMode);
                        }
                        else {
                            cerr << "Invalid recent files mode \""
                            << qPrintable(recentFilesModeName)
                            << "\" for \"-recent-files-mode\" option" << std::endl;
                            hasFatalError = true;
                        }
                    }
                    else {
                        cerr << "Missing recent files mode for \"-recent-files-mode\" option" << std::endl;
                        hasFatalError = true;
                    }
                } else if (thisParam == "-scene-load") {
                    if (myParams->hasNext()) {
                        myState.sceneFileName = myParams->nextString("Scene File Name");
                        if (myParams->hasNext()) {
                            myState.sceneNameOrNumber = myParams->nextString("Scene Name or Number");
                        }
                        else {
                            cerr << "Missing scene name/number for \"-scene-load\" option" << std::endl;
                            hasFatalError = true;
                        }
                    }
                    else {
                        cerr << "Missing scene file name for \"-scene-load\" option" << std::endl;
                        hasFatalError = true;
                    }
                } else if (thisParam == "-scene-load-hd") {
                    if (myParams->hasNext()) {
                        myState.sceneFileNameNoDialog = myParams->nextString("Scene File Name");
                        if (myParams->hasNext()) {
                            myState.sceneNameOrNumber = myParams->nextString("Scene Name or Number");
                        }
                        else {
                            cerr << "Missing scene name/number for " << thisParam << " option" << std::endl;
                            hasFatalError = true;
                        }
                    }
                    else {
                        cerr << "Missing scene file name for " << thisParam << " option" << std::endl;
                        hasFatalError = true;
                    }
                } else if (thisParam == "-spec-load-all") {
                    if ( ! myState.specFileNameLoadAll.isEmpty()) {
                        cerr << qPrintable(moreThanOneSpecFileErrorMessage) << endl;
                    }
                    else if (myParams->hasNext()) {
                        myState.specFileNameLoadAll = myParams->nextString("Spec File Name");
                    }
                    else {
                        cerr << "Missing spec file name for \"-spec\" option" << endl;
                        hasFatalError = true;
                    }
                } else if (thisParam == "-graphics-size") {
                    if (myParams->hasNext()) {
                        myState.graphicsSizeXY[0] = myParams->nextInt("Graphics Size X");
                    }
                    else {
                        cerr << "Missing X & Y sizes for graphics" << endl;
                        hasFatalError = true;
                    }
                    if (myParams->hasNext()) {
                        myState.graphicsSizeXY[1] = myParams->nextInt("Graphics Size Y");
                    }
                    else {
                        cerr << "Missing Y sizes for graphics" << endl;
                        hasFatalError = true;
                    }
                } else if (thisParam == "-window-size") {
                    if (myParams->hasNext()) {
                        myState.windowSizeXY[0] = myParams->nextInt("Window Size X");
                    }
                    else {
                        cerr << "Missing X & Y sizes for window" << endl;
                        hasFatalError = true;
                    }
                    if (myParams->hasNext()) {
                        myState.windowSizeXY[1] = myParams->nextInt("Window Size Y");
                    }
                    else {
                        cerr << "Missing Y sizes for window" << endl;
                        hasFatalError = true;
                    }
                } else if (thisParam == "-window-pos") {
                    if (myParams->hasNext()) {
                        myState.windowPosXY[0] = myParams->nextInt("Window Position X");
                    }
                    else {
                        cerr << "Missing X & Y position for window" << endl;
                        hasFatalError = true;
                    }
                    if (myParams->hasNext()) {
                        myState.windowPosXY[1] = myParams->nextInt("Window Position Y");
                    }
                    else {
                        cerr << "Missing Y position for window" << endl;
                        hasFatalError = true;
                    }
                } else if (thisParam.startsWith("-psn")) {
                    /*
                     * 21 April 2014 (Did not have this problem before this date)
                     *
                     * IGNORE this parameter.  For some reason, when a Mac
                     * version is started from Finder, a "-psn" parameter
                     * is being added to the parameters.  If this parameter
                     * is not ignored, Workbench starts, the icon bounces 
                     * a few times, and then Workbench quits (due to
                     * "unrecognized option", below), and the user is 
                     * not given any error message.
                     *
                     * http://stackoverflow.com/questions/10242115/os-x-strange-psn-command-line-parameter-when-launched-from-finder
                     * http://trac.wxwidgets.org/ticket/15432
                     */
                } else if (thisParam == "-NSDocumentRevisionsDebugMode") {
                    /*
                     * When wb_view is started within Apple's XCode, these parameters are added
                     * so ignore them:
                     *    -NSDocumentRevisionsDebugMode YES
                     */
                    if (myParams->hasNext()) {
                        myParams->nextString("Argument to -NSDocumentRevisionsDebugMode");
                    }
                } else {
                    cerr << "unrecognized option \"" << thisParam << "\"" << endl;
                    printHelp(progName);
                    hasFatalError = true;
                }
            } else {
                if (thisParam.endsWith(DataFileTypeEnum::toFileExtension(DataFileTypeEnum::SPECIFICATION))) {
                    if ( ! myState.specFileNameLoadWithDialog.isEmpty()) {
                        cerr << qPrintable(moreThanOneSpecFileErrorMessage) << endl;
                        hasFatalError = true;
                    }
                    else {
                        myState.specFileNameLoadWithDialog = thisParam;
                    }
                }
                else {
                    myState.fileList.push_back(thisParam);
                }
            }
        }
    }
    catch (const ProgramParametersException& e) {
        cerr << e.whatString() << std::endl;
        hasFatalError = true;
    }
    
    if ( ( ! myState.specFileNameLoadWithDialog.isEmpty())
        && ( ! myState.specFileNameLoadAll.isEmpty())) {
        cerr << qPrintable(moreThanOneSpecFileErrorMessage) << endl;
        hasFatalError = true;
    }
    
    if (hasFatalError) {
        exit(-1);
    }
    
    /*
     * If any files are listed, disable splash screen
     */
    if ( ! myState.fileList.empty()) {
        myState.showSplash = false;
    }
    if ( ! myState.sceneFileName.isEmpty()) {
        myState.showSplash = false;
    }
    if ( ! myState.sceneFileNameNoDialog.isEmpty()) {
        myState.showSplash = false;
    }
    if ( ! myState.specFileNameLoadWithDialog.isEmpty()) {
        myState.showSplash = false;
    }
    if ( ! myState.specFileNameLoadAll.isEmpty()) {
        myState.showSplash = false;
    }
}

ProgramState::ProgramState()
{
    sceneFileName = "";
    sceneFileNameNoDialog = "";
    sceneNameOrNumber = "";
    windowSizeXY[0] = -1;
    windowSizeXY[1] = -1;
    windowPosXY[0] = -1;
    windowPosXY[1] = -1;
    graphicsSizeXY[0] = -1;
    graphicsSizeXY[1] = -1;
    showSplash = true;
    macMenuFlag = false;
}
