#if defined(G__WIN32) && defined(__CINT__) && !defined(__MAKECINT__)
{
   Info("PrepareExportJson.C", "Has to be run in compiled mode ... doing this for you.");
   gSystem->CompileMacro("PrepareExportJson.C");
   PrepareExportJson();
}
#else

//#include "TUUID.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoBoolNode.h"
#include "TGeoMatrix.h"
#include "TGeoCompositeShape.h"
#include "TBuffer3D.h"
#include "TString.h"
#include "TColor.h"
#include "Riostream.h"
#include <algorithm>
#include <queue>
#include <math.h>
//TString doGuidStuff()
//{
//   TUUID uuid;
//   TString str_uuid = uuid.AsString();
//   str_uuid.ToUpper();
//   return StrDup(str_uuid);
//}


//====================GUID SYSTEM SELECTION BEGIN ======================
// PLEASE comment out the lines below to select your system code
//======================================================================
// this is the mac and ios version 
#define GUID_CFUUID
//======================================================================
// obviously this is the windows version
//#define GUID_WINDOWS
//======================================================================
// android version that uses a call to a java api
//#define GUID_ANDROID
//======================================================================
// This is the linux friendly implementation, but it could work on other
// systems that have libuuid available
//#define GUID_LIBUUID
//====================GUID SYSTEM SELECTION END ======================

//============== GUID =================begin
/*
The MIT License (MIT)

Copyright (c) 2014 Graeme Hill (http://graemehill.ca)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#include "guid.h"

#ifdef GUID_LIBUUID
#include <uuid/uuid.h>
#endif

#ifdef GUID_CFUUID
#include <CoreFoundation/CFUUID.h>
#endif

#ifdef GUID_WINDOWS
#include <objbase.h>
#endif

#ifdef GUID_ANDROID
#include <jni.h>
#endif

using namespace std;

// overload << so that it's easy to convert to a string
ostream &operator<<(ostream &s, const Guid &guid)
{
  return s << hex << setfill('0')
    << setw(2) << (int)guid._bytes[0]
    << setw(2) << (int)guid._bytes[1]
    << setw(2) << (int)guid._bytes[2]
    << setw(2) << (int)guid._bytes[3]
    << "-"
    << setw(2) << (int)guid._bytes[4]
    << setw(2) << (int)guid._bytes[5]
    << "-"
    << setw(2) << (int)guid._bytes[6]
    << setw(2) << (int)guid._bytes[7]
    << "-"
    << setw(2) << (int)guid._bytes[8]
    << setw(2) << (int)guid._bytes[9]
    << "-"
    << setw(2) << (int)guid._bytes[10]
    << setw(2) << (int)guid._bytes[11]
    << setw(2) << (int)guid._bytes[12]
    << setw(2) << (int)guid._bytes[13]
    << setw(2) << (int)guid._bytes[14]
    << setw(2) << (int)guid._bytes[15];
}

// create a guid from vector of bytes
Guid::Guid(const vector<unsigned char> &bytes)
{
  _bytes = bytes;
}

// create a guid from array of bytes
Guid::Guid(const unsigned char *bytes)
{
  _bytes.assign(bytes, bytes + 16);
}

// converts a single hex char to a number (0 - 15)
unsigned char hexDigitToChar(char ch)
{
  if (ch > 47 && ch < 58)
    return ch - 48;

  if (ch > 96 && ch < 103)
    return ch - 87;

  if (ch > 64 && ch < 71)
    return ch - 55;

  return 0;
}

// converts the two hexadecimal characters to an unsigned char (a byte)
unsigned char hexPairToChar(char a, char b)
{
  return hexDigitToChar(a) * 16 + hexDigitToChar(b);
}

// create a guid from string
Guid::Guid(const string &fromString)
{
  _bytes.clear();

  char charOne, charTwo;
  bool lookingForFirstChar = true;

  for (const char &ch : fromString)
  {
    if (ch == '-')
      continue;

    if (lookingForFirstChar)
    {
      charOne = ch;
      lookingForFirstChar = false;
    }
    else
    {
      charTwo = ch;
      auto byte = hexPairToChar(charOne, charTwo);
      _bytes.push_back(byte);
      lookingForFirstChar = true;
    }
  }

}

// create empty guid
Guid::Guid()
{
  _bytes = vector<unsigned char>(16, 0);
}

// copy constructor
Guid::Guid(const Guid &other)
{
  _bytes = other._bytes;
}

// overload assignment operator
Guid &Guid::operator=(const Guid &other)
{
  _bytes = other._bytes;
  return *this;
}

// overload equality operator
bool Guid::operator==(const Guid &other) const
{
  return _bytes == other._bytes;
}

// overload inequality operator
bool Guid::operator!=(const Guid &other) const
{
  return !((*this) == other);
}

// This is the linux friendly implementation, but it could work on other
// systems that have libuuid available
#ifdef GUID_LIBUUID
Guid GuidGenerator::newGuid()
{
  uuid_t id;
  uuid_generate(id);
  return id;
}
#endif

// this is the mac and ios version 
#ifdef GUID_CFUUID
Guid GuidGenerator::newGuid()
{
  auto newId = CFUUIDCreate(NULL);
  auto bytes = CFUUIDGetUUIDBytes(newId);
  CFRelease(newId);

  const unsigned char byteArray[16] =
  {
    bytes.byte0,
    bytes.byte1,
    bytes.byte2,
    bytes.byte3,
    bytes.byte4,
    bytes.byte5,
    bytes.byte6,
    bytes.byte7,
    bytes.byte8,
    bytes.byte9,
    bytes.byte10,
    bytes.byte11,
    bytes.byte12,
    bytes.byte13,
    bytes.byte14,
    bytes.byte15
  };
  return byteArray;
}
#endif

// obviously this is the windows version
#ifdef GUID_WINDOWS
Guid GuidGenerator::newGuid()
{
  GUID newId;
  CoCreateGuid(&newId);

  const unsigned char bytes[16] = 
  {
    (newId.Data1 >> 24) & 0xFF,
    (newId.Data1 >> 16) & 0xFF,
    (newId.Data1 >> 8) & 0xFF,
    (newId.Data1) & 0xff,

    (newId.Data2 >> 8) & 0xFF,
    (newId.Data2) & 0xff,

    (newId.Data3 >> 8) & 0xFF,
    (newId.Data3) & 0xFF,

    newId.Data4[0],
    newId.Data4[1],
    newId.Data4[2],
    newId.Data4[3],
    newId.Data4[4],
    newId.Data4[5],
    newId.Data4[6],
    newId.Data4[7]
  };

  return bytes;
}
#endif

// android version that uses a call to a java api
#ifdef GUID_ANDROID
GuidGenerator::GuidGenerator(JNIEnv *env)
{
  _env = env;
  _uuidClass = env->FindClass("java/util/UUID");
  _newGuidMethod = env->GetStaticMethodID(_uuidClass, "randomUUID", "()Ljava/util/UUID;");
  _mostSignificantBitsMethod = env->GetMethodID(_uuidClass, "getMostSignificantBits", "()J");
  _leastSignificantBitsMethod = env->GetMethodID(_uuidClass, "getLeastSignificantBits", "()J");
}

Guid GuidGenerator::newGuid()
{
  jobject javaUuid = _env->CallStaticObjectMethod(_uuidClass, _newGuidMethod);
  jlong mostSignificant = _env->CallLongMethod(javaUuid, _mostSignificantBitsMethod);
  jlong leastSignificant = _env->CallLongMethod(javaUuid, _leastSignificantBitsMethod);

  unsigned char bytes[16] = 
  {
    (mostSignificant >> 56) & 0xFF,
    (mostSignificant >> 48) & 0xFF,
    (mostSignificant >> 40) & 0xFF,
    (mostSignificant >> 32) & 0xFF,
    (mostSignificant >> 24) & 0xFF,
    (mostSignificant >> 16) & 0xFF,
    (mostSignificant >> 8) & 0xFF,
    (mostSignificant) & 0xFF,
    (leastSignificant >> 56) & 0xFF,
    (leastSignificant >> 48) & 0xFF,
    (leastSignificant >> 40) & 0xFF,
    (leastSignificant >> 32) & 0xFF,
    (leastSignificant >> 24) & 0xFF,
    (leastSignificant >> 16) & 0xFF,
    (leastSignificant >> 8) & 0xFF,
    (leastSignificant) & 0xFF,
  };
  return bytes;
}
#endif

TString UpperCase(TString str){
  for (int i=0;i<strlen(str);i++) str[i]=toupper(str[i]);
  return str;
}

TString doGuidStuff(GuidGenerator generator)
{
        auto myGuid = generator.newGuid();
        std::stringstream stream;
        stream << myGuid;
        auto guidString = stream.str();

        return UpperCase((TString)guidString);
}


//============== GUID ================= End


//============== JSON GEOMETRY ======== Begin

TString CheckVolume(TGeoShape *geoShape)
{
   // Check the type of volume and ...

   TString clsname = geoShape->ClassName();
   TString VolFound = clsname;

     bb=10;
   //process different shapes
   if (clsname == "TGeoBBox") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoParaboloid") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoSphere") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoConeSeg") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoCone") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoTubeSeg") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoTube") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoPcon") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoTorus") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoPgon") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoHype") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoScaledShape") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoArb8") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoPara") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoTrap") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoGtra") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoTrd1") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoTrd2") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoCtub") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoEltu") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoXtru") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
   } else if (clsname == "TGeoShapeAssembly") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
      //cout << "\nERROR! Skiping Solid NOT supported: "<<clsname<<"\n";
   } else if (clsname == "TGeoUnion") {
      //VolFound="";
      //cout << "\nERROR! Skiping Solid NOT supported: "<<clsname<<"\n";
   //} else if (clsname == "TGeoCompositeShape") {
      //cout << "\n---\n-> Geometry of volume: "<< clsname << " \n---\n";
      //cout << "\nERROR! Skiping Solid NOT supported: "<<clsname<<"\n";
   //} else if (clsname == "TGeoIntersection") {
      //VolFound="";
      //cout << "\nERROR! Skiping Solid NOT supported: "<<clsname<<"\n";
   //} else if (clsname == "TGeoSubtraction") {
      //   VolFound="";
      //cout << "\nERROR! Skiping Solid NOT supported: "<<clsname<<"\n";
   } else {
      VolFound = "";
      cout << "\nERROR! Skiping Solid NOT supported: " <<"\t" <<clsname<< "\t" <<geoShape->GetName()<< " \n";
   }
   return VolFound;
}



//______________________________________________________________________________

void CalculateSurfaceNormal(Double_t p1[3], Double_t p2[3], Double_t p3[3], Double_t *Vnormal)
{
//
// Here we compute two vectors coplanar to a face. The faces are triangles defined by points
// p1,p2,p3, then two possible vectors are: U = p2 - p1 and V = p3 - p1
// With the two vectors, U & V we compute the cross product between them to find a perpendicular
// vector to the face: Vn = U x V
// Vn(x) = Uy * Vz - Uz * Vy
// Vn(y) = Uz * Vx - Ux * Vz
// Vn(z) = Ux * Vy - Uy * Vx
//
//The length of the vector: Vmod=SQRT(Vnx**2+Vny**2+Vnz**2)
//
//
   Double_t U[3], V[3], Vn[3], Vmod;

   U[0] = p2[0] - p1[0];
   U[1] = p2[1] - p1[1];
   U[2] = p2[2] - p1[2];

   V[0] = p3[0] - p1[0];
   V[1] = p3[1] - p1[1];
   V[2] = p3[2] - p1[2];

   Vn[0] = U[1] * V[2] - U[2] * V[1];
   Vn[1] = U[2] * V[0] - U[0] * V[2];
   Vn[2] = U[0] * V[1] - U[1] * V[0];

   Vmod = sqrt(pow(Vn[0], 2.0) + pow(Vn[1], 2.0) + pow(Vn[2], 2.0));

   if (Vmod == 0) {
      //cout << "\n-->Error:====================================================================\n";
      //for (int i = 0; i < 3; i++) {
         //cout << "-->P: " << p1[i] << " " << p2[i] << " " << p3[i] << " U: " << U[i] << " V: " << V[i] << " Vn: " << Vn[i] << "\n";
      //}
      //cout << "\n-->==========================================================================\n";
      Vnormal[0] = 0;
      Vnormal[1] = 0;
      Vnormal[2] = 0;
   } else {
      Vnormal[0] = Vn[0] / Vmod;
      Vnormal[1] = Vn[1] / Vmod;
      Vnormal[2] = Vn[2] / Vmod;
      //cout << "Vmode:" << Vmod <<" Normal:" << Vnormal <<"\n";
   }
}


void CreateJsonVol(TGeoShape *shape, ofstream &JsonFile, TString guid)
{
   // Load some root geometry
   TBuffer3D *buffer_vol;
   int Nfaces = 0;
   queue <int *> Faces;
   queue <Double_t *> Fnormals;
   Double_t V1[3], V2[3], V3[3];
   Double_t *Vnormal;
   int *FaceTriangle;

   //--Consider each polygon of min 3 (triangles) and max 4 (quads) vertices.
   //--vrtx containes the index to vertexes for each segment of the polugon (max 2*4=8)
   //--vquad contains [continous] vertexes of the polygon (max=4)
   int vrtx[24], vquad[4];

   //     --Call the triangulation routine
   //cout << "\n-->Calling the triangulation routine---" << "\n";
   buffer_vol = shape->MakeBuffer3D();

   //cout << "\nMakeBuffer3D: No of Vert/Seg/poly: " << buffer_vol->NbPnts() << " / " << buffer_vol->NbSegs()<< " / " << buffer_vol->NbPols() <<" : "<<shape->GetName()<<"\n";

   //Start GEOMETRY definition in JSON with "{".
   JsonFile << "{\n";
   JsonFile << "\t\"uuid\": \"";
   //Generate a GUID for the specific volume for JSON file
   JsonFile << guid << "\",\n";
   JsonFile << "\t\"type\": \"Geometry\",\n";
   JsonFile << "\t\"data\":{\n";
   JsonFile << "\t\"vertices\": [";

   //Write all vertexes in JSON file
   int Kvertex = 0;
   for (int i = 0; i < buffer_vol->NbPnts() * 3; i = i + 3) {
      JsonFile << buffer_vol->fPnts[i] << ",  " << buffer_vol->fPnts[i + 1] << ",  " << buffer_vol->fPnts[i + 2] ;
      //cout << "--> " << i << " <-- " << buffer_vol->fPnts[i] << ",  " << buffer_vol->fPnts[i + 1] << ",  " << buffer_vol->fPnts[i + 2] << "\n";
      if (i < buffer_vol->NbPnts() * 3 - 3) {
         JsonFile << ","; //JsonFile <<",\n";
      } else {
         //JsonFile <<"\n";
      }
      Kvertex = Kvertex + 1;
   }
   JsonFile << "]," << "\n";
   //cout <<"Done vertex";


   Int_t IndexSegPolygon = 1; //define only triangles
   Nfaces = 0;
   //cout << "\n-->Generating triangles and normals ...";

   //      Loop through all the polygons that have been defined for this volume
   //      Each polygon will be converted to triangles faces saved in the que <Faces>
   //      and for each triangle we calculate the normal stored also in the queue <Fnormals>

   for (int Kpolygon = 0; Kpolygon < buffer_vol->NbPols(); Kpolygon++) {
      int iv = 0;

      //          Loop through all the segments of the current polygon and get the Vertex indexes
      //          expected max 4 segments and consequently 8 vertex indexes
      //          fSegs:  c0, p0, q0, c1, p1, q1, ..... ..... ....
      //          fPols;  c0, n0, s0, s1, ... sn, c1, n1, s0, ... sn
      if (buffer_vol->fPols[IndexSegPolygon]<=4){
        for (int k = 1; k < buffer_vol->fPols[IndexSegPolygon] + 1; k++) {
          vrtx[iv] = buffer_vol->fSegs[3 * buffer_vol->fPols[IndexSegPolygon + k] + 1];
          //cout <<"["<< vrtx[iv] <<",";
          iv = iv + 1;
          vrtx[iv] = buffer_vol->fSegs[3 * buffer_vol->fPols[IndexSegPolygon + k] + 2];
          //cout << vrtx[iv] <<"]";
          iv = iv + 1;
        }
        //cout << "\n";
        //          Segments are defined in the correct order for the definition of the normals
        //          but for each segment the sequence of vertexes may be inverted.
        //          check continuity of vertexes in the consequtive segments
        //          If Seg1 made of P0-P1 and Seg2 made of P0-P2 then Seg1 should be defined as P1-P0

        if (vrtx[1] != vrtx[2]) {
          if (vrtx[0] == vrtx[2]) {
            int itemp = vrtx[0];
            vrtx[0] = vrtx[1];
            vrtx[1] = itemp;

          } else if (vrtx[0] == vrtx[3]) {
            int itemp = vrtx[0];
            vrtx[0] = vrtx[1];
            vrtx[1] = itemp;

            itemp = vrtx[2];
            vrtx[2] = vrtx[3];
            vrtx[3] = itemp;

          } else if (vrtx[1] == vrtx[3]) {
            int itemp = vrtx[2];
            vrtx[2] = vrtx[3];
            vrtx[3] = itemp;
          }
        }

        //          ---
        //          Out of the (max 8) vertexes definining the (max 4) segments,
        //          get the (max 4) vertexes defining the polygon (in quad)
        //          Check that segments and vertexes are in the correct order
        //          Starting from the first segment (vrtx[0]-vrtx[1]),
        //          the second segment must have the vertex vrtx[1] and this should be the first in the order
        //          ---

        vquad[0] = vrtx[0];
        vquad[1] = vrtx[1];
        vrtx[0] = -1;
        vrtx[1] = -1;

        int kvert = 2;
        while (kvert < buffer_vol->fPols[IndexSegPolygon]) {
          for (int ks = 1; ks < iv; ks++) {
            if (vrtx[ks * 2 + 1] == vquad[kvert - 1]) {
               vquad[kvert] = vrtx[ks * 2];
               vrtx[ks * 2] = -1;
               vrtx[ks * 2 + 1] = -1;
               kvert = kvert + 1;
               break;
            } else if (vrtx[ks * 2] == vquad[kvert - 1]) {
               vquad[kvert] = vrtx[ks * 2 + 1];
               vrtx[ks * 2] = -1;
               vrtx[ks * 2 + 1] = -1;
               kvert = kvert + 1;
               break;
            }
          }
        }


        //     Using the vertex indexes get the x,y,z coordinates of each vertex of the 1st triangle

        for (int i = 0; i < 3; i++) {
          V1[i] = buffer_vol->fPnts[3 * vquad[2] + i];
          V2[i] = buffer_vol->fPnts[3 * vquad[1] + i];
          V3[i] = buffer_vol->fPnts[3 * vquad[0] + i];
        }

        if ((V1[0] == V2[0] && V1[1] == V2[1] && V1[2] == V2[2]) ||
            (V1[0] == V3[0] && V1[1] == V3[1] && V1[2] == V3[2]) ||
            (V2[0] == V3[0] && V2[1] == V3[1] && V2[2] == V3[2])) {

          //cout << "\nSkiping triangle face1: <" << Nfaces << "> [" << V1[0] << "," << V1[1] << " ," << V1[2] << "]";

        } else {
          //            Get the first triangle of the defined polygon
          FaceTriangle = new int[3];
          FaceTriangle[0] = vquad[2];
          FaceTriangle[1] = vquad[1];
          FaceTriangle[2] = vquad[0];

          Vnormal = new Double_t[3];
          CalculateSurfaceNormal(V1, V2, V3, Vnormal);

          Fnormals.push(Vnormal);

          Faces.push(FaceTriangle);
          Nfaces = Nfaces + 1;
        }

        //          If the polygon is a quad, get also the second triangle of the defined polygon
        //
        //         3 _____ 2
        //          |    /|
        //          | B / |
        //          |  /  |
        //          | / A |
        //          |/____|
        //         0       1
        //
        //         First trianle A: 0,1,2 second triangle B: 2,3,0

        if (buffer_vol->fPols[IndexSegPolygon] == 4) {
          //            Using the vertex indexes get the x,yz coordinates of each vertex of the 2nd triangle
          for (int i = 0; i < 3; i++) {
            V1[i] = buffer_vol->fPnts[3 * vquad[0] + i];
            V2[i] = buffer_vol->fPnts[3 * vquad[3] + i];
            V3[i] = buffer_vol->fPnts[3 * vquad[2] + i];
          }

          if ((V1[0] == V2[0] && V1[1] == V2[1] && V1[2] == V2[2]) ||
              (V1[0] == V3[0] && V1[1] == V3[1] && V1[2] == V3[2]) ||
              (V2[0] == V3[0] && V2[1] == V3[1] && V2[2] == V3[2])) {
            //cout << "\nSkiping triangle face2: <" << Nfaces << "> [" << V1[0] << "," << V1[1] << " ," << V1[2] << "]";

          } else {

            FaceTriangle = new int[3];
            FaceTriangle[0] = vquad[0];
            FaceTriangle[1] = vquad[3];
            FaceTriangle[2] = vquad[2];

            Vnormal = new Double_t[3];
            CalculateSurfaceNormal(V1, V2, V3, Vnormal);

            Fnormals.push(Vnormal);

            Faces.push(FaceTriangle);

            Nfaces = Nfaces + 1;
          }
        }
      }
      IndexSegPolygon = IndexSegPolygon + buffer_vol->fPols[IndexSegPolygon] + 2;
   } // end loop through polygons

   //      print the Normals (one normal per triangle)

   JsonFile << "\t\"normals\": [";
   //cout << "\n-->Generated "<< Nfaces << " Normals ... Writing normals to file ... \n";
   int Knormal = 0;
   while (!Fnormals.empty()) {
      Double_t *Vnorm = Fnormals.front();
      //cout << "\nNormal coordinates" << Vnorm[0] << "," << Vnorm[1] << "," << Vnorm[2];
      JsonFile << Vnorm[0] << "," << Vnorm[1] << "," << Vnorm[2];
      Knormal = Knormal + 1;

      if (Nfaces != Knormal) {
         JsonFile << ",";
      }
      //cout << "delete Normal pointer "<<Knormal<< " / "<<Nfaces<<" \n";

      Fnormals.pop();
      delete [] Vnorm;
   }
   JsonFile << " ],\n";

   // Start writing faces definition in Json file

   JsonFile << "\t\"faces\": [";

   //cout << "Done Normals\n";
   //cout << "-->Creating faces ... \n";

   //      print the triangles / faces definition
   int Kfaces = 0;
   while (!Faces.empty()) {
      int *point = Faces.front();
      JsonFile << "16," << point[0] << "," << point[1] << "," << point[2] << "," << Kfaces << "  ";

      Kfaces = Kfaces + 1;

      if (Nfaces != Kfaces) {
         JsonFile << ",";
      }
      //cout << "delete face pointer "<<Kfaces<< " / "<<Nfaces<<" \n";

      Faces.pop();
      delete [] point;
   }
   //cout << "Done Faces\n";
   JsonFile << "]";
   //Terminate the DATA definition for this volume with "}".
   JsonFile << "}\n";

   //Terminate the GEOMETRY definition for this volume with "}".
   JsonFile << "}";
}

TGeoShape *GetLshape(TGeoCompositeShape *shape)
{
   return (shape->GetBoolNode()->GetLeftShape());
}

TGeoShape *GetRshape(TGeoCompositeShape *shape)
{
   return (shape->GetBoolNode()->GetRightShape());
}


int GetJsonVol(TGeoNode *current,TGeoShape *shape, ofstream &JsonFile,TGeoShape *MapVol[],TString  *MapGuidVol,TString *MapGuidMat,int ColourIndex,int &nPhyVol, int kVolume){
  bool notfound;
  if (shape->IsComposite()) {
         kVolume=GetJsonVol(current,GetLshape((TGeoCompositeShape *)shape),JsonFile,MapVol,MapGuidVol,MapGuidMat,ColourIndex,nPhyVol,kVolume);
         kVolume=GetJsonVol(current,GetRshape((TGeoCompositeShape *)shape),JsonFile,MapVol,MapGuidVol,MapGuidMat,ColourIndex,nPhyVol,kVolume);
         nPhyVol++;
  } else if (shape->IsAssembly()) {
  } else {
    notfound=true;
    for (int i = 0; i < kVolume; i++) {
      if (shape == MapVol[i]) {
        notfound=false;
        break;
      }
    }
    if (notfound) {
      MapVol[kVolume] = shape;
      MapGuidVol[kVolume] = doGuidStuff(GuidGenerator());;

      MapGuidMat[current->GetColour()]="Used color";
      //if (kVolume > 0) 
      JsonFile << ",\n";
      CreateJsonVol(shape, JsonFile, MapGuidVol[kVolume]);
      kVolume++;
    }
  }
  return kVolume;
}




//============== JSON Materials definition ======== Begin



void ExportMaterials(int ColorIndex, TString guid, ofstream &JsonFile)
{
  Float_t r,g,b;
  int MatColor;

  TColor *color = gROOT->GetColor(ColorIndex);
  TString MatName=color->GetName();
  TString ColorHex=color->AsHexString();
  ColorHex.Remove(0,1);
  //MatColor=ColorHex.Atoi();

  //cout << "\nColor:" << ColorIndex<<"\t"<<ColorHex<<"\t";

  MatColor =  std::stoul((const std::string)ColorHex, nullptr, 16);

  //cout << MatColor <<"\n";

 
   JsonFile << "{\n\t\"uuid\": \"" << guid << "\",\n";
   JsonFile << "\t\"type\": \"MeshLambertMaterial\",\n";
   JsonFile << "\t\"name\":\"" << MatName << "\",\n";

   JsonFile << "\t\"color\":"<<MatColor<<",\n";
   JsonFile << "\t\"emissive\":"<<MatColor<<",\n";

   //if (MatName == "Vacuum" || (nA == 0 && nZ == 0)) {
   //   JsonFile << "\t\"opacity\": 0.0,\n";
   //   JsonFile << "\t\"transparent\": true\n }";
   //} else {
      JsonFile << "\t\"opacity\": 0.7,\n";
      JsonFile << "\t\"transparent\": true\n}";
   //}
}

//============== JSON Materials definition ======== End

void TabLine(int ntab, ofstream &JsonFile)
{
   for (int i = 0; i < ntab; i++) {
      JsonFile << "\t";
   }
}


//============== JSON Nodes definition ======== Begin



void SetTabInJson(int CurrentLevel,int PreviousLevel,ofstream &JsonFile){


  if (CurrentLevel > PreviousLevel) {
      //   If the level is greater that current level, then we define a new child
      JsonFile << ",\n";
      TabLine(CurrentLevel, JsonFile);
      JsonFile << "\"children\": [\n";
      TabLine(CurrentLevel, JsonFile);
      JsonFile << "{\n";
      //JsonFile<<"\n"<<"//       ("<<CurrentLevel<<","<<PreviousLevel<<")\n";
  } else if (CurrentLevel < PreviousLevel) {
      int nLev = PreviousLevel - CurrentLevel;
      JsonFile << "\n";
      for (int i = 0; i < nLev; i++) {
         TabLine(PreviousLevel - i, JsonFile);
         JsonFile << "}]\n";
      }
      TabLine(CurrentLevel, JsonFile);
      JsonFile << "},\n";
      TabLine(CurrentLevel, JsonFile);
      JsonFile<< "{\n";
      //JsonFile <<"//       ("<<CurrentLevel<<","<<PreviousLevel<<")\n";
  } else {
      JsonFile << "\n";
      TabLine(CurrentLevel, JsonFile);
      JsonFile << "},\n";
      TabLine(CurrentLevel, JsonFile);
      JsonFile << "{\n";
      //JsonFile <<"//       ("<<CurrentLevel<<","<<PreviousLevel<<")\n";
  }
}






TString GetHexRGB(int r, int g, int b)
{
   char hexcol[16];

   snprintf(hexcol, sizeof hexcol, "%02x%02x%02x", r, g, b);
   return hexcol;
}

TGeoCompositeShape * Conv(TGeoCompositeShape *shape){
  return shape;
}

Int_t ExportCurrentVolume(TGeoNode *current,TGeoShape *shape,TString guidMat,TGeoShape *MapVol[],TString *MapGuidVol,Int_t kVolume, Int_t nVol,  const Double_t * MatrixRot,const Double_t * MatrixTrans, int CurrentLevel, int PreviousLevel,TString Visible, ofstream &JsonFile)
{
  TGeoCompositeShape *Cshape;
  TString guid, guidVol;
  Int_t Kvol=0;
 
  if (shape->IsComposite()) {
    //JsonFile << "\n//          ....composite\n";
    Cshape=Conv((TGeoCompositeShape *)shape);
    MatrixRot   = Cshape->GetBoolNode()->GetLeftMatrix()->GetRotationMatrix() ;
    MatrixTrans = Cshape->GetBoolNode()->GetLeftMatrix()->GetTranslation();
    nVol=ExportCurrentVolume(current,Cshape->GetBoolNode()->GetLeftShape(),guidMat,MapVol,MapGuidVol,kVolume,nVol, MatrixRot, MatrixTrans, CurrentLevel, PreviousLevel,Visible, JsonFile);
    JsonFile << "\n";
    TabLine(CurrentLevel, JsonFile);
    JsonFile << "},\n";
    TabLine(CurrentLevel, JsonFile);
    JsonFile << "{\n";
    MatrixRot   = Cshape->GetBoolNode()->GetRightMatrix()->GetRotationMatrix() ;
    MatrixTrans = Cshape->GetBoolNode()->GetRightMatrix()->GetTranslation();
    nVol=ExportCurrentVolume(current,Cshape->GetBoolNode()->GetRightShape(),guidMat,MapVol,MapGuidVol,kVolume,nVol, MatrixRot, MatrixTrans, CurrentLevel, PreviousLevel,Visible, JsonFile);
  } else {
    nVol++;
    //   find the guid of the current volume
    if (shape->IsAssembly()) {
      guidVol = MapGuidVol[0];
    } else {
      guidVol = "";
      //cout << nVol<<"Shape:"<<shape->ClassName()<<"\t"<<shape->GetName()<<"\t";
      for (int i = 0; i < kVolume; i++) {
        if (shape == MapVol[i]) {
          guidVol = MapGuidVol[i];

          Kvol=i;
          break;
        }
      }
    }
    if (guidVol == "") {
      cout << "\n****ERROR: :[ "<<kVolume<<"\t"<< nVol<<"\t"<<guidVol<<"\tShape:"<<shape<<" ["<<current->GetNumber()<<"] "<<shape->ClassName()<<"\t"<<shape->GetName()<<"\t" << "  Not found **** STOP\n";
      return nVol;
    } 
    guid = doGuidStuff(GuidGenerator());

    TabLine(CurrentLevel + 1, JsonFile);
    JsonFile << "\"uuid\": \"" << guid << "\",\n";
    TabLine(CurrentLevel + 1, JsonFile);
    JsonFile << "\"type\": \"Mesh\",\n";
    TabLine(CurrentLevel + 1, JsonFile);
    JsonFile << "\"name\": \"" << current->GetVolume()->GetName()<<"-"<< current->GetNumber()<< "\","<<"\n";
    TabLine(CurrentLevel + 1, JsonFile);
    JsonFile << "\"geometry\": \"" << guidVol << "\",\n";
    TabLine(CurrentLevel + 1, JsonFile);
    JsonFile << "\"material\": \"" << guidMat << "\",\n";
    TabLine(CurrentLevel + 1, JsonFile);
//    JsonFile << "\"visible\": " << Visible << ",\n";
    JsonFile << "\"visible\": true,\n";
    TabLine(CurrentLevel + 1, JsonFile);
    JsonFile << "\"matrix\": [";

    for (int i = 0; i < 3; i++) {
     for (int l = 0; l < 3; l++) {
        JsonFile << MatrixRot[3 * l + i] << ",";
     }
     JsonFile << "0,";
    }

    for (int i = 0; i < 3; i++) {
      JsonFile << MatrixTrans[i] << ",";
    }
    JsonFile << "1 ]";
  }
  return nVol;
}



//============== Export JSON Tree ======== Begin



int CheckJsonVol(TGeoNode *current,TGeoShape *shape, ofstream &JsonFile,TGeoShape *MapVol[],TString  *MapGuidVol,TString *MapGuidMat,int ColourIndex, int &nPhyVol,int kVolume){
bool notfound;
   if (shape->IsComposite()) {
         kVolume=CheckJsonVol(current,GetLshape((TGeoCompositeShape *)shape),JsonFile,MapVol,MapGuidVol,MapGuidMat,ColourIndex,nPhyVol,kVolume);
         kVolume=CheckJsonVol(current,GetRshape((TGeoCompositeShape *)shape),JsonFile,MapVol,MapGuidVol,MapGuidMat,ColourIndex,nPhyVol,kVolume);
   } else {
    nPhyVol++;
    notfound=true;
    for (int i = 0; i < kVolume; i++) {
      if (shape == MapVol[i]) {
        notfound=false;
        break;
      }
    }
    if (notfound) {
      MapVol[kVolume] = shape;
      MapGuidVol[kVolume] = doGuidStuff(GuidGenerator());;

      MapGuidMat[current->GetColour()]="Used color";

      //JsonFile << kVolume << "\t"<<MapVol[kVolume]<< "\t"<<shape->ClassName()<<"\t"<< shape->GetName()<<"\t"<<current->GetNumber()<<"\t"<<current->GetColour()<<"\n";

      kVolume++;
    }
  }
  return kVolume;
}



int CheckLogicalVolumes(TGeoVolume *top, int MaxVisiLevel, int nAsblVol[], int StartingLevel,int &PreviousLevel, int kVolume,int nnd,TString *MapGuidMat,TGeoShape *MapVol[],TString *MapGuidVol,int &nPhyVol,bool OneCopy,ofstream &JsonFile){

  TGeoNode *current;
  TGeoVolume *vol;
  TGeoShape *shape;
  int kLevel = 0;
  int vstep=nnd/50;
  int nSkipTab1=0;
  int CurrentLevel=0;
  int RelLevel=0;
  int nIterVol=0;

 

  JsonFile << "Pre\tSta\tCur\tSkp\tLev\tRLV\tCName\tname\tCopy\n";


  TGeoIterator iter(top);
  int ncount=0;
  while ((current = iter.Next())) {
    ncount++;
    RelLevel=StartingLevel+iter.GetLevel();
    nAsblVol[RelLevel]=0;

    vol = current->GetVolume();
    shape = vol->GetShape();

    nSkipTab1=0;
    for (int j = 0; j<RelLevel; j++) {nSkipTab1+=nAsblVol[j];}
    CurrentLevel=RelLevel+nSkipTab1;

    JsonFile <<PreviousLevel<<"\t"<<StartingLevel<<"\t"<<CurrentLevel<<"\t"<<nSkipTab1<<"\t"<<iter.GetLevel()<<"\t"<<RelLevel;
    for (int j = 0; j<RelLevel; j++) {JsonFile << "\t-->";}
    JsonFile <<shape->ClassName()<<"\t"<<current->GetVolume()->GetName()<<"\t"<<current->GetNumber()<<"\n";

    if (vol->IsAssembly()) {
      nAsblVol[RelLevel]=-1;
    }
      if (CurrentLevel<=MaxVisiLevel){

      TString chvol = CheckVolume(shape);
      if (chvol != "") {
        bool notfound=true;
        for (int i = 0; i < kVolume; i++) {
          if (shape == MapVol[i]) {
            notfound=false;
            break;
          }
        }
        if (notfound) {
          kVolume=CheckJsonVol(current,shape,JsonFile,MapVol,MapGuidVol,MapGuidMat,current->GetColour(),nPhyVol,kVolume);
        }
        PreviousLevel=CurrentLevel;
      }
    }
    if (OneCopy&&(current->GetNumber()>1)) {iter.Skip();}
  }
  return(kVolume);
}


int ExportLogicalVolumes(TGeoVolume *top, int MaxVisiLevel, int nAsblVol[], int StartingLevel,int PreviousLevel, int kVolume,int nnd,TString *MapGuidMat,TGeoShape *MapVol[],TString *MapGuidVol,int &nPhyVol,bool OneCopy, ofstream &JsonFile){

  TGeoNode *current;
  TGeoVolume *vol;
  TGeoShape *shape;
  int kLevel = 0;
  int vstep=nnd/50;
  int nSkipTab1=0;
  int CurrentLevel=0;
  int RelLevel=0;
  int nIterVol=0;

   
  TGeoIterator iter(top);
  while ((current = iter.Next())) {
    RelLevel=StartingLevel+iter.GetLevel();
    nAsblVol[RelLevel]=1;

    vol = current->GetVolume();
    shape = vol->GetShape();

    nSkipTab1=0;
    for (int j = 0; j<RelLevel; j++) {nSkipTab1+=nAsblVol[j];}
    CurrentLevel=RelLevel+nSkipTab1;
    
    if (vol->IsAssembly()) {
      nAsblVol[RelLevel]=-1;
    }
    if (CurrentLevel<=MaxVisiLevel){

      TString chvol = CheckVolume(shape);
      if (chvol != "") {
        nPhyVol++;
        bool notfound=true;
        for (int i = 0; i < kVolume; i++) {
          if (shape == MapVol[i]) {
            notfound=false;
            break;
          }
        }
        if (notfound) {
          kVolume=GetJsonVol(current,shape,JsonFile,MapVol,MapGuidVol,MapGuidMat,current->GetColour(),nPhyVol,kVolume);
        }
        PreviousLevel=CurrentLevel;
        if (OneCopy&&(current->GetNumber()>1)) {iter.Skip();}
      } else {
          iter.Skip();
      }
    }
  }
  return(kVolume);
}



int ExportPhysicalVolumes(TGeoVolume *top, int MaxVisiLevel, int nAsblVol[], int StartingLevel,int &PreviousLevel, int nn,int nnd,TString * MapGuidMat,TGeoShape *MapVol[],TString *MapGuidVol,int kVolume, bool OneCopy,TString Vguid,ofstream &JsonFile){

  TGeoNode *current;
  TGeoVolume *vol;
  TGeoShape *shape;
  int vstep=nnd/50;
  int nSkipTab1=0;
  int CurrentLevel=0;
  int RelLevel=0;
  TString guidMat;
  TString MatName;
  TString Visible;
  
  if (top->IsAssembly()) {
  //  logfile << "ReStart: " <<"\t"<< StartingLevel<<"\t"<<PreviousLevel<<"\t"<<top->GetName()<<"\n";
  } else {
  //  logfile << "Start: " <<"\t Normal volume"<<"\t"<<top->GetName()<< "\t"<<top->GetNumber()<<"\n";
  }

  TGeoIterator iter(top);
  int ncount=0;
  while ((current = iter.Next())) {
    ncount++;
    RelLevel=StartingLevel+iter.GetLevel();
    nAsblVol[RelLevel]=0;

    vol = current->GetVolume();
    shape = vol->GetShape();

    nSkipTab1=0;
    for (int j = 0; j<RelLevel; j++) {nSkipTab1+=nAsblVol[j];}

    CurrentLevel=RelLevel+nSkipTab1;

    if (vol->IsAssembly()) {
      nAsblVol[RelLevel]=-1;
    }
    if (CurrentLevel<=MaxVisiLevel){

    TString chvol = CheckVolume(shape);
    if (chvol != "") {
      //   Get the guid of the material of the current volume
      guidMat = MapGuidMat[current->GetColour()];
         
      SetTabInJson(CurrentLevel,PreviousLevel,JsonFile);
      //Visible="false";
      //if (vol->IsVisible()){Visible="true";}
      //if (vol->IsAssembly()|| !vol->IsVisible()){guidMat=Vguid;}
      MatName=vol->GetMaterial()->GetName();
      if (vol->IsAssembly()||MatName=="dummy"){guidMat=Vguid;}

      if (current->GetMotherVolume()->IsAssembly()){
        const Double_t *MatrixRot   = current->GetMatrix()->GetRotationMatrix();
        const Double_t *MatrixTrans = current->GetMatrix()->GetTranslation();
        nn=ExportCurrentVolume(current,shape,guidMat, MapVol,MapGuidVol,kVolume, nn,  MatrixRot, MatrixTrans, CurrentLevel, PreviousLevel,Visible, JsonFile);
      } else {
        const Double_t *MatrixRot   = current->GetMatrix()->GetRotationMatrix();
        const Double_t *MatrixTrans = current->GetMatrix()->GetTranslation();
        nn=ExportCurrentVolume(current,shape,guidMat, MapVol,MapGuidVol,kVolume, nn,  MatrixRot, MatrixTrans, CurrentLevel, PreviousLevel,Visible, JsonFile);
      }
          
      PreviousLevel=CurrentLevel;
      if (OneCopy&&(current->GetNumber()>1)) {iter.Skip();}
    } else {
      iter.Skip();
    }
    
  }
  JsonFile << "\n";

  return(nn);
}


int ExportJsonTree(TGeoVolume *top, int MaxVisiLevel, int nAsblVol[], int StartingLevel,int &PreviousLevel, int kVolume,int nnd,TString *MapGuidMat,TGeoShape *MapVol[],TString *MapGuidVol,int &nPhyVol,bool OneCopy,ofstream &JsonTree){

  TGeoNode *current;
  TGeoVolume *vol;
  TGeoShape *shape;
  int kLevel = 0;
  int vstep=nnd/50;
  int nSkipTab1=0;
  int CurrentLevel=0;
  int RelLevel=0;
  int nIterVol=0;

  TGeoIterator iter(top);
  int ncount=0;
  while ((current = iter.Next())) {
    ncount++;
    RelLevel=StartingLevel+iter.GetLevel();
    nAsblVol[RelLevel]=0;

    vol = current->GetVolume();
    shape = vol->GetShape();

    nSkipTab1=0;
    for (int j = 0; j<RelLevel; j++) {nSkipTab1+=nAsblVol[j];}
    CurrentLevel=RelLevel+nSkipTab1;


    if (CurrentLevel<=MaxVisiLevel){
      JsonTree <<iter.GetLevel()<<"\t";
      for (int j = 0; j<RelLevel; j++) {JsonTree << "-->";}
      JsonTree <<current->GetVolume()->GetName()<<"\t"<<current->GetNumber()<<"\t"<<shape->ClassName()<<"\t"<<vol->GetMaterial()->GetName()<<"\n";
    }
     if (vol->IsAssembly()) {
      nAsblVol[RelLevel]=-1;
    }
    if (CurrentLevel<=MaxVisiLevel){

      TString chvol = CheckVolume(shape);
      if (chvol != "") {
        bool notfound=true;
        for (int i = 0; i < kVolume; i++) {
          if (shape == MapVol[i]) {
            notfound=false;
            break;
          }
        }
        if (notfound) {
          //if (current->IsVisible()) {cout <<"\n NOT VISIBLE " << vol->GetName() <<"\t"<<shape->ClassName()<<"\t"<<vol->GetMaterial()->GetName()<<"\t"<< current->IsVisible() <<"\n";}

          kVolume=CheckJsonVol(current,shape,JsonTree,MapVol,MapGuidVol,MapGuidMat,current->GetColour(),nPhyVol,kVolume);
        }
        PreviousLevel=CurrentLevel;

      } else {
        iter.Skip();
      }
        
    }
    if (OneCopy&&(current->GetNumber()>1)) {iter.Skip();}
  }
  return(kVolume);
}

//============== Export JSON Tree ======== End




//============== JSON Nodes definition ======== End

void PrepareExportJson()
{
   cout << "\n JSON export functions loaded ....\n";
}

#endif
