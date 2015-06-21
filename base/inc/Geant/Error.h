//===------------------ Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantError.h
 * @brief Error handling routines.
 * @author Philippe Canal
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_ERROR_H
#define GEANT_ERROR_H

#ifndef GEANT_CONFIG_H
#include "Geant/Config.h"
#endif

namespace Geant {
   enum class EMsgLevel {
      kUnset    =  -1,
      kPrint    =   0,
      kInfo     =   1000,
      kWarning  =   2000,
      kError    =   3000,
      kBreak    =   4000,
      kSysError =   5000,
      kFatal    =   6000
   };

   inline namespace GEANT_IMPL_NAMESPACE {

#ifndef GEANT_NVCC
      void ErrorHandlerImpl(EMsgLevel level, const char *location, const char *msgfmt, ...);
#endif

      template <typename... ArgsTypes>
      GEANT_CUDA_BOTH_CODE
      void MessageHandler(EMsgLevel level, const char *location, const char *msgfmt, ArgsTypes... params) {
#ifdef GEANT_NVCC
         const char *type = nullptr;
         switch(level) {
            case EMsgLevel::kPrint: type = "Print"; break;
            case EMsgLevel::kInfo: type = "Info"; break;
            case EMsgLevel::kWarning: type = "Warning"; break;
            case EMsgLevel::kError: type = "Error"; break;
            case EMsgLevel::kBreak: type = "Break"; break;
            case EMsgLevel::kSysError: type = "SysError"; break;
            case EMsgLevel::kFatal: type = "Fatal"; break;
            default: type = "Unknown Level"; break;
         }
         if (level == EMsgLevel::kPrint)
            printf("%s:",location);
         else
            printf("%s in <%s>:",type, location);
         printf(msgfmt,params...);
         printf("\n");
         if (level >= EMsgLevel::kFatal) {
#ifdef GEANT_CUDA_DEVICE_BUILD
            // Did not find a way to halt a kernel from within yet.
            //cudaDeviceReset();
            //cudaThreadExit();
            //throw("Fatal error in CUDA kernel");
#else
            exit( EXIT_FAILURE );
#endif
         }
#else
         // Currently we use the ROOT message handler on the host/gcc code.
         Geant::cxx::ErrorHandlerImpl(level,location,msgfmt,params...);
#endif
      }

      template <typename... ArgsTypes>
      GEANT_CUDA_BOTH_CODE
      void Warning(const char *location, const char *msgfmt, ArgsTypes... params)
      {
         MessageHandler(EMsgLevel::kWarning,location,msgfmt, params...);
      }

      template <typename... ArgsTypes>
      GEANT_CUDA_BOTH_CODE
      void Error(const char *location, const char *msgfmt, ArgsTypes... params)
      {
         MessageHandler(EMsgLevel::kError,location,msgfmt, params...);
      }

      template <typename... ArgsTypes>
      GEANT_CUDA_BOTH_CODE
      void Fatal(const char *location, const char *msgfmt, ArgsTypes... params)
      {
         MessageHandler(EMsgLevel::kFatal,location,msgfmt, params...);
      }

   } // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_ERROR_H
