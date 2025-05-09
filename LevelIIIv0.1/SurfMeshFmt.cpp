
#ifdef MEMDEBUG
#include "MemDbg.hpp"
#endif

#include "SurfMeshFmt.hpp"
#include "FileIterator.hpp"
#include "Set.hpp"

#include <cstdio>
#include <iomanip>

using std::setw ;
using std::setprecision ;
using std::endl ;
using std::ios ;
using std::istream ;
using std::ostream ;

/* ------------------------------------------------------------------------

   Define the SurfMeshReader, which is a helper class

*/

class SurfMeshReader {

    public:

        SurfMeshReader(
            FileIterator& iin,
            FemModel& model) ;

        void ReadFile() ;

    private:

        FemModel& model ;
        FileIterator& in ;

        Dict<String,void (SurfMeshReader::*)()> dispatch ;

        Dict<int,int> FacetRgn0 ;
        Dict<int,int> FacetRgn1 ;
        Dict<int,int> RegionMat ;
        Dict<int,int> Retained ;

        void Vertex() ;
        void Facet() ;
        void Region() ;

} ;

/* --------------------------------------------------------------------

    These are the reader and writer functions

*/

bool ReadSurfMeshFmt(
    istream& in,
    FemModel& model)
{
    FileIterator iter(in) ;
    SurfMeshReader reader(iter,model) ;
    reader.ReadFile() ;
    return(true) ;
}

/* ---------------------------------------------------------------------

   Here define the methods for the SurfMeshReader class

*/

SurfMeshReader::SurfMeshReader(
    FileIterator& iin,
    FemModel &imodel) :
        model(imodel),
        in(iin)
{    
    dispatch.Store(String("VTX:"),    &SurfMeshReader::Vertex) ;
    dispatch.Store(String("FACET:"),  &SurfMeshReader::Facet) ;
    dispatch.Store(String("REGION:"), &SurfMeshReader::Region) ;
}

void SurfMeshReader::ReadFile()
{
    // parse data and dispatch to appropriate routines

    for (in.First() ; in.More() ; ++in) {
        if ((in.Len() > 0) && (dispatch.HasKey(in[0]))) {
            void (SurfMeshReader::*handler)() = dispatch[in[0]] ;
            (this->*handler)() ;
        }
    }

    // create model attributes facet regions and region to material mapping

    model.ModelAtts.AddIntMapAtt("FacetRgn0",FacetRgn0) ;
    model.ModelAtts.AddIntMapAtt("FacetRgn1",FacetRgn1) ;
    model.ModelAtts.AddIntMapAtt("RegionMat",RegionMat) ;
    model.ModelAtts.AddIntMapAtt("RetainedFacets",Retained) ;
}

void SurfMeshReader::Vertex()
{
    // this handles node definitions:
    //   VTX: id x y z

    model.AddNode(in[1].AsInt(),
                  Vec3D(in[2].AsFloat(),
                        in[3].AsFloat(),
                        in[4].AsFloat()),
                  Vec3D(0,0,0)) ;
}

void SurfMeshReader::Facet()
{
    // this handles an facet definitions:
    //   FACET: id num_nodes n0 n1 ... region0 region1 retained_flag

    int eid = in[1].AsInt() ;
    int nnodes = in[2].AsInt() ;

    List<int> conn ;
    for (int i=0 ; i<nnodes ; ++i) conn.Append(in[i+3].AsInt()) ;

    FemModel::ElemType e_type ;
    int rgn0, rgn1 ;
    int retained ;

    if (nnodes == 4) {
        e_type = FemModel::SURF_QUAD_4 ;
        rgn0 = in[7].AsInt() ;
        rgn1 = in[8].AsInt() ;
        retained = (in[9].AsInt() != 0) ;
    } else {
        e_type = FemModel::SURF_TRI_3 ;
        rgn0 = in[6].AsInt() ;
        rgn1 = in[7].AsInt() ;
        retained = (in[8].AsInt() != 0) ;
    }

    model.AddElem(eid,1,0,e_type,conn) ;
    FacetRgn0.Store(eid,rgn0) ;
    FacetRgn1.Store(eid,rgn1) ;
    Retained.Store(eid,retained) ;
}

void SurfMeshReader::Region()
{
    // this handles a region to material mapping definitions:
    //   REGION: region material

    RegionMat.Store(in[1].AsInt(),in[2].AsInt()) ;
}


