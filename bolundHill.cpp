/* The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <cmath>
#include <memory>
#include "complexDynamics/dynamicSmagorinskyDynamics.h"

using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
//#define DESCRIPTOR descriptors::MRTD3Q19Descriptor
#define DESCRIPTOR descriptors::D3Q19Descriptor

plint extraLayer      = 0;  // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
const plint blockSize = 20; // Zero means: no sparse representation.
const plint envelopeWidth = 1;  // For standard BGK dynamics.
const plint extendedEnvelopeWidth = 2;  // Because the Guo off lattice boundary condition
                                        //   needs 2-cell neighbor access.

bool performOutput = false;
bool doImages = false;
bool useAllDirections = false;
bool useRegularizedWall = false;
bool useIncompressible = false;
bool poiseuilleInlet = false;
bool convectiveScaling = false;

T kinematicViscosity       = 0.;
T averageInletVelocity     = 0.;
plint referenceResolution  = 0.;
T nuLB                     = 0.;
T fluidDensity             = 0.;
T volume                   = 0.;
T userDefinedInletDiameter = 0.;

plint referenceDirection = 0;
plint openingSortDirection = 0;

T simTime = 0;
plint startLevel = 0;
plint maxLevel   = 0;
T epsilon = 0.;

TriangleSet<T>* triangleSet = 0;
T currentTime = 0;

class InletProfile {
public: 
    InletProfile(T dx_, T dt_) : dx(dx_), dt(dt_)
    { }
    void operator() (plint iX, plint iY, plint iZ, Array<T,3>& u) const {
        T u_str = 0.7866 * dt / dx;
        T rndCoeff = ((double)rand() / (double)RAND_MAX) - 0.5;
        T z_0 = 3.e-4 / dx;
        T u_0 = 0.4 * dt / dx;
        T kappa = 0.4;
        T z = (T)iZ / dx;
        // u[0]  = u_0 / kappa * log(z/z_0);
        u[0]  = u_0 / kappa * log(z/z_0) + u_str*rndCoeff;
        //u[0]  = u_0 / kappa * log((z+1.e-9)/z_0);
        u[1]  = T();
        u[2]  = T();
    }
private:
    T dx;
    T dt;
};


// Structure which defines an ``opening''. The surface geometry of the aneurysm,
//   as given by the user in the form of an STL file, contains holes, which in 
//   the specific simulation represent inlets and outlets.
template<typename T>
struct Opening {
    bool inlet;
    Array<T,3> center;
    T innerRadius;
};

std::vector<Opening<T> > openings;

void iniLattice( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 VoxelizedDomain3D<T>& voxelizedDomain, T uAveLB)
{
    // Switch all remaining outer cells to no-dynamics, except the outer
    //   boundary layer, and keep the rest as BGKdynamics.
    defineDynamics(lattice, voxelizedDomain.getVoxelMatrix(), lattice.getBoundingBox(),
                   new NoDynamics<T,DESCRIPTOR>, voxelFlag::outside);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T) 1., Array<T,3>(uAveLB,(T) 0.,(T) 0.));
    lattice.initialize();
}

// This function assigns proper boundary conditions to the openings of the surface geometry
//   of the aneurysm. Which opening is inlet and which is outlet is defined by the user in
//   the input XML file. For the inlet, there is a choice between a Poiseuille velocity
//   profile and a simple plug velocity profile. At the outlets a Neumann boundary
//   condition with constant pressure is prescribed.
void setOpenings (
    std::vector<BoundaryProfile3D<T,Velocity>*>& inletOutlets,
    TriangleBoundary3D<T>& boundary, T uLB, T dx, T dt )
{
    for (pluint i=0; i<openings.size(); ++i) {
        Opening<T>& opening = openings[i];
        opening.center = computeBaryCenter (
                boundary.getMesh(),
                boundary.getInletOutlet(openingSortDirection)[i] );
        opening.innerRadius = computeInnerRadius (
                boundary.getMesh(),
                boundary.getInletOutlet(openingSortDirection)[i] );

        if (opening.inlet) {
            if (poiseuilleInlet) {
                inletOutlets.push_back (
                        new PoiseuilleProfile3D<T>(uLB) );
                    //new WindProfile3D<T>(uLB, dx, dt) ); ////wahahahahagaggasgduasgdgas
            }
            else {
                inletOutlets.push_back (
                        new VelocityPlugProfile3D<T>(uLB) );
            }
        }
        else {
            inletOutlets.push_back (
                    new DensityNeumannBoundaryProfile3D<T> );
        }
    }
}

// This function outputs velocity, vorticity and pressure data, at selected
//   points of the computational domain, given their coordinates in physical units.
void winddata( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
         Array<T,3> location, T i, T dx, T dt ) 
{
    plb_ofstream measurePoints("velocities.dat", std::ofstream::app);
    std::vector< Array<T,3> > physicalPositions, positions;

    physicalPositions.push_back(Array<T,3>(-82.8,  29.5000,  6.0500)); //M0Z05S
    physicalPositions.push_back(Array<T,3>(-82.8,  28.7000, 13.0500)); //M0Z12S

    physicalPositions.push_back(Array<T,3>( 45.6, 101.0,  2.9000)); //M1Z02S
    physicalPositions.push_back(Array<T,3>( 45.6, 101.0,  5.9000)); //M1Z05S
    physicalPositions.push_back(Array<T,3>( 45.6, 101.0,  9.8000)); //M1Z09S

    physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 11.9000)); //M2Z01S
    physicalPositions.push_back(Array<T,3>( 63.1000, 111.8000, 12.9000)); //M2Z02S
    physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 14.4000)); //M2Z03S
    physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 15.9000)); //M2Z05S
    physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 19.9000)); //M2Z09S

    physicalPositions.push_back(Array<T,3>(101.2, 132.0, 13.7000)); //M3Z02S
    physicalPositions.push_back(Array<T,3>(101.2, 132.0, 16.7000)); //M3Z05S
    physicalPositions.push_back(Array<T,3>(101.2, 132.0, 20.7000)); //M3Z09S

    physicalPositions.push_back(Array<T,3>(149.5, 162.6,  2.8000)); //M4Z02S
    physicalPositions.push_back(Array<T,3>(149.5, 162.6,  5.8000)); //M4Z05S
    physicalPositions.push_back(Array<T,3>(149.5, 162.6,  9.8000)); //M4Z09S

    physicalPositions.push_back(Array<T,3>( 99.5,  83.1,  4.8000)); //M5Z02S
    physicalPositions.push_back(Array<T,3>( 99.5,  83.1,  7.8000)); //M5Z05S

    physicalPositions.push_back(Array<T,3>( 51.9, 132.2, 13.4000)); //M6Z02S
    physicalPositions.push_back(Array<T,3>( 51.9, 132.2, 16.4000)); //M6Z05S

    physicalPositions.push_back(Array<T,3>( 31.1, 132.0,  2.8000)); //M7Z02S
    physicalPositions.push_back(Array<T,3>( 31.1, 132.0,  5.8000)); //M7Z05S

    physicalPositions.push_back(Array<T,3>(190.0, 131.9,  3.8000)); //M8Z02S
    physicalPositions.push_back(Array<T,3>(190.0, 131.9,  6.7000)); //M8Z05S

    physicalPositions.push_back(Array<T,3>(425.3,  93.6000,  6.4000)); //M9Z05S

    physicalPositions.push_back(Array<T,3>(-82.8,  29.3000,  3.0500)); //M0Z02C
    physicalPositions.push_back(Array<T,3>(-82.8,  28.7,  6.0500)); //M0Z05C
    physicalPositions.push_back(Array<T,3>(-82.8,  28.7,  10.0500)); //M0Z09C
    physicalPositions.push_back(Array<T,3>(-82.8,  30.3000,  16.0500)); //M0Z15C

    physicalPositions.push_back(Array<T,3>( 63.2, 110.8900, 19.9000)); //M2Z09C
    physicalPositions.push_back(Array<T,3>( 63.1000, 112.7000, 21.9000)); //M2Z11C

    physicalPositions.push_back(Array<T,3>(101.1000, 135.6000, 20.7000)); //M3Z09C

    physicalPositions.push_back(Array<T,3>( 51.9, 132.2000, 20.4000)); //M6Z09C

    physicalPositions.push_back(Array<T,3>(190.0, 131.9000,  10.8000)); //M8Z09C

    physicalPositions.push_back(Array<T,3>(425.3,  92.7,  3.3000)); //M9Z02C
    physicalPositions.push_back(Array<T,3>(425.3,  92.7,  6.4000)); //M9Z05C
    physicalPositions.push_back(Array<T,3>(425.3,  92.7,  10.4000)); //M9Z09C
    physicalPositions.push_back(Array<T,3>(425.3,  92.7,  17.)); //M9Z15C


    positions.resize(physicalPositions.size());
    T currentSimTime = i*dt;
    measurePoints << i << "; " << currentSimTime ;
    for (pluint k=0; k<physicalPositions.size(); ++k) {
        positions[k] = (physicalPositions[k]-location)/dx;
    }

    std::vector<Array<T,3> > velocities = velocitySingleProbes(lattice, positions);
    std::vector<Array<T,3> > vorticities = vorticitySingleProbes(lattice, positions);
    std::vector<T> densities = densitySingleProbes(lattice, positions);

    std::vector<T> data;
    for (pluint k=0; k<physicalPositions.size(); ++k) {
        Array<T,3> pos = physicalPositions[k];
        Array<T,3> vel = velocities[k]*dx/dt;
        Array<T,3> vort = vorticities[k]/dt;
        
        T pressure = DESCRIPTOR<T>::cs2*(densities[k]-1.)*dx*dx/(dt*dt)*fluidDensity;
//         if (performOutput) {
//             pcout << "Pos ("
//                   << pos[0] << "," << pos[1] << "," << pos[2]
//                   << "); Velocity ("
//                   << vel[0] << "," << vel[1] << "," << vel[2]
//                   << "); Vorticity ("
//                   << vort[0] << "," << vort[1] << "," << vort[2]
//                   << "); Pressure " << pressure << std::endl;
//         }
        if (performOutput) {
               measurePoints 
                     << "; " << pos[0] << "; " << pos[1] << "; " << pos[2] 
                     << "; " << vel[0] << "; " << vel[1] << "; " << vel[2]
                     << "; " << vort[0] << "; " << vort[1] << "; " << vort[2];
        }
            
        data.push_back(norm(vel));
        data.push_back(norm(vort));
        data.push_back(pressure);
    }
    
    measurePoints << std::endl;
}

void windprofile( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
         Array<T,3> location, T i, T dx, T dt ) 
{
    T xOffset = 82.8;

    plb_ofstream measureProfiles0("profile-0.dat", std::ofstream::app);
    plb_ofstream measureProfiles1("profile-1.dat", std::ofstream::app);
    plb_ofstream measureProfiles2("profile-2.dat", std::ofstream::app);
    plb_ofstream measureProfiles3("profile-3.dat", std::ofstream::app);
    plb_ofstream measureProfiles4("profile-4.dat", std::ofstream::app);
    plb_ofstream measureProfiles5("profile-5.dat", std::ofstream::app);
    plb_ofstream measureProfiles6("profile-6.dat", std::ofstream::app);
    plb_ofstream measureProfiles7("profile-7.dat", std::ofstream::app);
    plb_ofstream measureProfiles8("profile-8.dat", std::ofstream::app);
    plb_ofstream measureProfiles9("profile-9.dat", std::ofstream::app);
    std::vector<Array<T,3> > physicalPositions, positions;
    
    physicalPositions.push_back(Array<T,3>(-82.8+xOffset,  28.7000,  6.0500)); //M0
    physicalPositions.push_back(Array<T,3>( 45.6+xOffset,    101.0,  6.0500)); //M1
    physicalPositions.push_back(Array<T,3>( 63.2+xOffset,    110.9,  6.0500)); //M2
    physicalPositions.push_back(Array<T,3>(101.2+xOffset,    132.0,  6.0500)); //M3
    physicalPositions.push_back(Array<T,3>(149.5+xOffset,    162.6,  6.0500)); //M4
    physicalPositions.push_back(Array<T,3>( 99.5+xOffset,     83.1,  6.0500)); //M5
    physicalPositions.push_back(Array<T,3>( 51.9+xOffset,    132.2,  6.0500)); //M6
    physicalPositions.push_back(Array<T,3>( 31.1+xOffset,    132.0,  6.0500)); //M7
    physicalPositions.push_back(Array<T,3>(190.0+xOffset,    131.9,  6.0500)); //M8
    physicalPositions.push_back(Array<T,3>(425.3+xOffset,     92.7,  6.0500)); //M9
    // physicalPositions.push_back(Array<T,3>(-82.8,  29.5000,  6.0500)); //M0Z05S
    // physicalPositions.push_back(Array<T,3>(-82.8,  28.7000, 13.0500)); //M0Z12S

    // physicalPositions.push_back(Array<T,3>( 45.6, 101.0,  2.9000)); //M1Z02S
    // physicalPositions.push_back(Array<T,3>( 45.6, 101.0,  5.9000)); //M1Z05S
    // physicalPositions.push_back(Array<T,3>( 45.6, 101.0,  9.8000)); //M1Z09S

    // physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 11.9000)); //M2Z01S
    // physicalPositions.push_back(Array<T,3>( 63.1000, 111.8000, 12.9000)); //M2Z02S
    // physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 14.4000)); //M2Z03S
    // physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 15.9000)); //M2Z05S
    // physicalPositions.push_back(Array<T,3>( 63.2, 110.9, 19.9000)); //M2Z09S

    // physicalPositions.push_back(Array<T,3>(101.2, 132.0, 13.7000)); //M3Z02S
    // physicalPositions.push_back(Array<T,3>(101.2, 132.0, 16.7000)); //M3Z05S
    // physicalPositions.push_back(Array<T,3>(101.2, 132.0, 20.7000)); //M3Z09S

    // physicalPositions.push_back(Array<T,3>(149.5, 162.6,  2.8000)); //M4Z02S
    // physicalPositions.push_back(Array<T,3>(149.5, 162.6,  5.8000)); //M4Z05S
    // physicalPositions.push_back(Array<T,3>(149.5, 162.6,  9.8000)); //M4Z09S

    // physicalPositions.push_back(Array<T,3>( 99.5,  83.1,  4.8000)); //M5Z02S
    // physicalPositions.push_back(Array<T,3>( 99.5,  83.1,  7.8000)); //M5Z05S

    // physicalPositions.push_back(Array<T,3>( 51.9, 132.2, 13.4000)); //M6Z02S
    // physicalPositions.push_back(Array<T,3>( 51.9, 132.2, 16.4000)); //M6Z05S

    // physicalPositions.push_back(Array<T,3>( 31.1, 132.0,  2.8000)); //M7Z02S
    // physicalPositions.push_back(Array<T,3>( 31.1, 132.0,  5.8000)); //M7Z05S

    // physicalPositions.push_back(Array<T,3>(190.0, 131.9,  3.8000)); //M8Z02S
    // physicalPositions.push_back(Array<T,3>(190.0, 131.9,  6.7000)); //M8Z05S

    // physicalPositions.push_back(Array<T,3>(425.3,  93.6000,  6.4000)); //M9Z05S

    // physicalPositions.push_back(Array<T,3>(-82.8,  29.3000,  3.0500)); //M0Z02C
    // physicalPositions.push_back(Array<T,3>(-82.8,  28.7,  6.0500)); //M0Z05C
    // physicalPositions.push_back(Array<T,3>(-82.8,  28.7,  10.0500)); //M0Z09C
    // physicalPositions.push_back(Array<T,3>(-82.8,  30.3000,  16.0500)); //M0Z15C

    // physicalPositions.push_back(Array<T,3>( 63.2, 110.8900, 19.9000)); //M2Z09C
    // physicalPositions.push_back(Array<T,3>( 63.1000, 112.7000, 21.9000)); //M2Z11C

    // physicalPositions.push_back(Array<T,3>(101.1000, 135.6000, 20.7000)); //M3Z09C

    // physicalPositions.push_back(Array<T,3>( 51.9, 132.2000, 20.4000)); //M6Z09C

    // physicalPositions.push_back(Array<T,3>(190.0, 131.9000,  10.8000)); //M8Z09C

    // physicalPositions.push_back(Array<T,3>(425.3,  92.7,  3.3000)); //M9Z02C
    // physicalPositions.push_back(Array<T,3>(425.3,  92.7,  6.4000)); //M9Z05C
    // physicalPositions.push_back(Array<T,3>(425.3,  92.7,  10.4000)); //M9Z09C
    // physicalPositions.push_back(Array<T,3>(425.3,  92.7,  17.)); //M9Z15C


    positions.resize(physicalPositions.size());
    T currentSimTime = i*dt;
    // measureProfiles << i << "; " << currentSimTime ;
    for (pluint k=0; k<physicalPositions.size(); ++k) {
        positions[k] = (physicalPositions[k]-location)/dx;
    }

    // std::vector<T> profiles = *computeVelocityNorm(lattice, line0)
    // plint nx = lattice.getNx();
    // plint ny = lattice.getNy();
    plint nz = lattice.getNz();
    // Box3D line0(100, 100, 100, 100, 0,                   nz);
    // Box3D line1(50, 50, 100, 100, 0,                   nz);

    Box3D line0( positions[0][0], positions[0][0], positions[0][1], positions[0][1], 0, nz);
    Box3D line1( positions[1][0], positions[1][0], positions[1][1], positions[1][1], 0, nz);
    Box3D line2( positions[2][0], positions[2][0], positions[2][1], positions[2][1], 0, nz);
    Box3D line3( positions[3][0], positions[3][0], positions[3][1], positions[3][1], 0, nz);
    Box3D line4( positions[4][0], positions[4][0], positions[4][1], positions[4][1], 0, nz);
    Box3D line5( positions[5][0], positions[5][0], positions[5][1], positions[5][1], 0, nz);
    Box3D line6( positions[6][0], positions[6][0], positions[6][1], positions[6][1], 0, nz);
    Box3D line7( positions[7][0], positions[7][0], positions[7][1], positions[7][1], 0, nz);
    Box3D line8( positions[8][0], positions[8][0], positions[8][1], positions[8][1], 0, nz);
    Box3D line9( positions[9][0], positions[9][0], positions[9][1], positions[9][1], 0, nz);
    
    if (performOutput) {

        measureProfiles0 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line0) 
            << *computeVelocityComponent(lattice, line0, 0) 
            << *computeVelocityComponent(lattice, line0, 1) 
            << *computeVelocityComponent(lattice, line0, 2) << std::endl;

            measureProfiles1 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line2) 
            << *computeVelocityComponent(lattice, line1, 0) 
            << *computeVelocityComponent(lattice, line1, 1) 
            << *computeVelocityComponent(lattice, line1, 2) << std::endl;

            measureProfiles2 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line2) 
            << *computeVelocityComponent(lattice, line2, 0) 
            << *computeVelocityComponent(lattice, line2, 1) 
            << *computeVelocityComponent(lattice, line2, 2) << std::endl;
			
            measureProfiles3 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line3) 
            << *computeVelocityComponent(lattice, line3, 0) 
            << *computeVelocityComponent(lattice, line3, 1) 
            << *computeVelocityComponent(lattice, line3, 2) << std::endl;

            measureProfiles4 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line4) 
            << *computeVelocityComponent(lattice, line4, 0) 
            << *computeVelocityComponent(lattice, line4, 1) 
            << *computeVelocityComponent(lattice, line4, 2) << std::endl;

            measureProfiles5 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line5) 
            << *computeVelocityComponent(lattice, line5, 0) 
            << *computeVelocityComponent(lattice, line5, 1) 
            << *computeVelocityComponent(lattice, line5, 2) << std::endl;

            measureProfiles6 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line6) 
            << *computeVelocityComponent(lattice, line6, 0) 
            << *computeVelocityComponent(lattice, line6, 1) 
            << *computeVelocityComponent(lattice, line6, 2) << std::endl;

            measureProfiles7 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line7) 
            << *computeVelocityComponent(lattice, line7, 0) 
            << *computeVelocityComponent(lattice, line7, 1) 
            << *computeVelocityComponent(lattice, line7, 2) << std::endl;

            measureProfiles8 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line8) 
            << *computeVelocityComponent(lattice, line8, 0) 
            << *computeVelocityComponent(lattice, line8, 1) 
            << *computeVelocityComponent(lattice, line8, 2) << std::endl;

            measureProfiles9 << i << " " << currentSimTime << " " << dx << " " << dt << " " << nz << " " 
            << *computeVelocityNorm(lattice, line0) 
            << *computeVelocityComponent(lattice, line9, 0) 
            << *computeVelocityComponent(lattice, line9, 1) 
            << *computeVelocityComponent(lattice, line9, 2) << std::endl;

            
        // measureProfiles << positions << "; " ;
    }
    
    // measureProfiles << std::endl;
}

void writeVTK (
         MultiBlockLattice3D<T,DESCRIPTOR>& lattice, 
         MultiScalarField3D<int>& flagMatrix,
         std::string fname, Array<T,3> location, T dx, T dt )
{  
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();
    plint nz = lattice.getNz();


    // Box3D a_plane(0, nx-1, util::roundToInt(132./dx-1.), util::roundToInt(132./dx+1.), 0, nz-1);
    Box3D b_plane(0, nx-1, util::roundToInt(132./dx-1.), util::roundToInt(132./dx+1.), 0, nz-1);

    VtkImageOutput3D<T> vtkOut(fname, dx, location);
///     Physical Units
//     vtkOut.writeData<float>(*computeDensity(lattice, vtkDomain), "p", util::sqr(dx/dt)*fluidDensity);
//     vtkOut.writeData<float>(*computeVelocityNorm(lattice, vtkDomain), "u", dx/dt);
//     vtkOut.writeData<float>(*copyConvert<int,T>(*extractSubDomain(flagMatrix, vtkDomain)), "flag", 1.);

///     Lattice Untits    
    // vtkOut.writeData<float>(*computeDensity(lattice, vtkDomain), "p", 1.);

    vtkOut.writeData<float>(*computeVelocityNorm(lattice, b_plane), "velocityNorm", dx/dt);
    vtkOut.writeData<3, float>(*computeVelocity(lattice, b_plane), "velocity", dx/dt);
    // vtkOut.writeData<float>(*computeVelocityComponent(lattice, b_plane, 0), "u", dx/dt);
    // vtkOut.writeData<float>(*computeVelocityComponent(lattice, b_plane, 1), "v", dx/dt);
    // vtkOut.writeData<float>(*computeVelocityComponent(lattice, b_plane, 2), "w", dx/dt);
    vtkOut.writeData<float>(*copyConvert<int,T>(*extractSubDomain(flagMatrix, b_plane)), "flag", 1.);
    
}

void writeFullVTK (
         MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
         MultiScalarField3D<int>& flagMatrix,
         Box3D const& vtkDomain, std::string fname, Array<T,3> location, T dx, T dt )
{  

    VtkImageOutput3D<T> vtkOut(fname, dx, location);
///     Physical Units
//     vtkOut.writeData<float>(*computeDensity(lattice, vtkDomain), "p", util::sqr(dx/dt)*fluidDensity);
//     vtkOut.writeData<float>(*computeVelocityNorm(lattice, vtkDomain), "u", dx/dt);
//     vtkOut.writeData<float>(*copyConvert<int,T>(*extractSubDomain(flagMatrix, vtkDomain)), "flag", 1.);

///     Lattice Untits    
    vtkOut.writeData<float>(*computeDensity(lattice, vtkDomain), "p", 1.);
    // vtkOut.writeData<float>(*computeVelocityNorm(lattice, vtkDomain), "u", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice, vtkDomain), "velocity", dx/dt);
    vtkOut.writeData<float>(*copyConvert<int,T>(*extractSubDomain(flagMatrix, vtkDomain)), "flag", 1.);
    
}

void writeImages (
         MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
         MultiScalarField3D<int>& flagMatrix,
         Box3D const& vtkDomain, std::string fname, Array<T,3> location, T dx, T dt )
{
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();
    plint nz = lattice.getNz();
    
    Box3D x_plane(util::roundToInt(nx/2), util::roundToInt(nx/2), 0, ny-1, 0, nz-1);
    Box3D y_plane(0, nx-1, util::roundToInt(ny/2), util::roundToInt(ny/2), 0, nz-1);
    Box3D z_plane(0, nx-1, 0, ny-1, util::roundToInt(nz/2), util::roundToInt(nz/2));
    
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledPpm(fname+'_'+util::val2str(0), *computeVelocityNorm(lattice, x_plane));
    imageWriter.writeScaledPpm(fname+'_'+util::val2str(1), *computeVelocityNorm(lattice, y_plane));
    imageWriter.writeScaledPpm(fname+'_'+util::val2str(2), *computeVelocityNorm(lattice, z_plane));    
}

// void createSpongeZone (MultiBlockLattice3D<T,DESCRIPTOR>& lattice, T uLB, T dx, T dt)
// {
//     T bulkValue;
//     plint numSpongeCells;
//     T outletSpongeZoneWidth;
//     Array<plint,6> numSPongeCells;
//     
//     
//     bulkValue = lattice.getOmega;
//     outletSpongeZoneWidth = 0.4;
//     numOutletSpongeCells = util::roundToInt(outletSpongeZoneWidth/dx);
//     
//     
//     
//         // Number of sponge zone lattice nodes at all the outer domain boundaries.
//         // So: 0 means the boundary at x = 0
//         //     1 means the boundary at x = nx-1
//         //     2 means the boundary at y = 0
//         //     and so on...
//         numSpongeCells[0] = 0;
//         numSpongeCells[1] = param.numOutletSpongeCells;
//         numSpongeCells[2] = 0;
//         numSpongeCells[3] = 0;
//         numSpongeCells[4] = 0;
//         numSpongeCells[5] = 0;
//     
//     std::vector<MultiBlock3D*> args;
//     args.push_back(lattice);
//     applyProcessingFunctional(new ViscositySpongeZone3D<T,DESCRIPTOR>(
//         lattice.getNx, lattice.getNy, lattice.getNz, bulkValue, numSpongeCells),
//         lattice->getBoundingBox(),args);
// }

// This is the function that prepares and performs the actual simulation.
std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > run (
        plint level, MultiBlockLattice3D<T,DESCRIPTOR>* iniVal=0 )
{
    plint margin = 3; // Extra margin of allocated cells around the obstacle. 
    plint borderWidth = 1; // Because the Guo boundary condition acts in a one-cell layer.
                           // Requirement: margin>=borderWidth.

    // The resolution is doubled at each coordinate direction with the increase of the
    //   resolution level by one. The parameter ``referenceResolution'' is by definition
    //   the resolution at grid refinement level 0.
    plint resolution = referenceResolution * util::twoToThePower(level);

    // The next few lines of code are typical. They transform the surface geometry of the
    //   aneurysm given by the user to more efficient data structures that are internally
    //   used by palabos. The TriangleBoundary3D structure will be later used to assign
    //   proper boundary conditions.
    DEFscaledMesh<T>* defMesh =
        new DEFscaledMesh<T>(*triangleSet, resolution, referenceDirection, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    // When convective scaling is used (relationship of dt with respect to dx as the grid is
    //   refined) the value of the kinematic viscosity must be also properly adjusted.
    T nuLB_ = nuLB;
    if (convectiveScaling) {
        nuLB_ = nuLB * util::twoToThePower(level);
    }
    T dx = boundary.getDx();
    T dt = nuLB_ / kinematicViscosity *dx*dx;
    T uAveLB = averageInletVelocity *dt/dx;
    T omega = 1./(3.*nuLB_+0.5);
    Array<T,3> location(boundary.getPhysicalLocation());

    std::string output = "myOutput";
    pcout << "uLB=" << uAveLB << std::endl;
    pcout << "nuLB=" << nuLB_ << std::endl;
    pcout << "tau=" << 1./omega << std::endl;
    if (performOutput) {
        pcout << "dx=" << dx << std::endl;
        pcout << "dt=" << dt << std::endl;
    }

    // The aneurysm simulation is an interior (as opposed to exterior) flow problem. For
    //   this reason, the lattice nodes that lay inside the computational domain must
    //   be identified and distinguished from the ones that lay outside of it. This is
    //   handled by the following voxelization process.
    const int flowType = voxelFlag::inside;
    VoxelizedDomain3D<T> voxelizedDomain (
            boundary, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize );
    if (performOutput) {
        pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    }

    MultiScalarField3D<int> flagMatrix((MultiBlock3D&)voxelizedDomain.getVoxelMatrix());
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::inside, flagMatrix.getBoundingBox(), 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::innerBorder, flagMatrix.getBoundingBox(), 1);
    pcout << "Number of fluid cells: " << computeSum(flagMatrix) << std::endl;

    Dynamics<T,DESCRIPTOR>* dynamics = 0;
    if (useIncompressible) {
        dynamics = new IncBGKdynamics<T,DESCRIPTOR>(omega); // In this model velocity equals momentum.
    }
    else {
        //dynamics = new BGKdynamics<T,DESCRIPTOR>(omega); // In this model velocity equals momentum
                                                         //   divided by density.
		//dynamics = new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(omega, 0.14);
        dynamics = new ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T,DESCRIPTOR>(omega, 0.14);
    }
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice 
        = generateMultiBlockLattice<T,DESCRIPTOR> (
                voxelizedDomain.getVoxelMatrix(), envelopeWidth, dynamics );
    lattice->toggleInternalStatistics(true);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    
    // Array<T,3> inletRealPos(-82.79, 125., 20.); // inlet slightly displaced in the pipe
    // Array<T,3> outlet1RealPos(439.99, 125., 20.); // outlet slightly displaced in the pipe
    Array<T,3> inletRealPos(0.1, 125., 20.); // inlet slightly displaced in the pipe
    Array<T,3> outlet1RealPos(522, 125., 20.); // outlet slightly displaced in the pipe
    
    //Array<T,3> outlet2RealPos(0.007,0.06,0.105);
    //T diameterReal = 0.0612;
//    T diameterReal = 0.0622;
    T diameterReal = 30.;
    T domainHeightReal = 40.;
    T domainWidthReal  = 250.;
    
    
    Array<T,3> inletPos((inletRealPos-location)/dx);
    Array<T,3> outlet1Pos((outlet1RealPos-location)/dx);
    //Array<T,3> outlet2Pos((outlet2RealPos-location)/dx);
    plint diameter = util::roundToInt(diameterReal/dx);
    plint domainHeight= util::roundToInt(domainHeightReal/dx);
    plint domainWidth = util::roundToInt(domainWidthReal/dx);
    
    Box3D inletDomain(util::roundToInt(inletPos[0]), util::roundToInt(inletPos[0]),
                      util::roundToInt(inletPos[1]-domainWidth), util::roundToInt(inletPos[1]+domainWidth),
                      util::roundToInt(inletPos[2]-domainHeight), util::roundToInt(inletPos[2]+domainHeight));
    Box3D behindInlet(inletDomain.x0-diameter, inletDomain.x0-1,
                      inletDomain.y0, inletDomain.y1,
                      inletDomain.z0, inletDomain.z1);

    Box3D outlet1Domain(util::roundToInt(outlet1Pos[0]), util::roundToInt(outlet1Pos[0]),
                        util::roundToInt(outlet1Pos[1]-domainWidth), util::roundToInt(outlet1Pos[1]+domainWidth),
                        util::roundToInt(outlet1Pos[2]-domainHeight), util::roundToInt(outlet1Pos[2]+domainHeight));
    Box3D behindOutlet1(outlet1Domain.x0+1, outlet1Domain.x0+diameter,
                        outlet1Domain.y0, outlet1Domain.y1,
                        outlet1Domain.z0, outlet1Domain.z1);

    boundaryCondition->addVelocityBoundary0N(inletDomain, *lattice);
        
    setBoundaryVelocity(*lattice, inletDomain, 
        InletProfile(dx, dt));
        //Array<T,3>(uAveLB,(T)0.,(T)0.) );
//     boundaryCondition->addPressureBoundary0P(outlet1Domain, *lattice);
//     setBoundaryDensity(*lattice, outlet1Domain, (T)1.);
    
 
    
//    boundaryCondition->addPressureBoundary0N(outlet2Domain, *lattice);
//    setBoundaryDensity(*lattice, outlet2Domain, (T)1.);

    defineDynamics(*lattice, flagMatrix, lattice->getBoundingBox(), new BounceBack<T,DESCRIPTOR>(1.), 0);
    defineDynamics(*lattice, behindInlet, new BounceBack<T,DESCRIPTOR>(1.));
    defineDynamics(*lattice, behindOutlet1, new BounceBack<T,DESCRIPTOR>(1.));
//    defineDynamics(*lattice, behindOutlet2, new BounceBack<T,DESCRIPTOR>(1.));
    
       //free slip top
    plint nx = lattice->getNx();
    plint ny = lattice->getNy();
    plint nz = lattice->getNz();
    
    Box3D water(1, nx-2, 1, ny-2, 3, 3);
    Box3D top(1, nx-2, 1, ny-2, nz-4, nz-4);
    Box3D rightWall(1, nx-2, 3, 3, 0, nz-1);
    Box3D leftWall(1, nx-2, ny-4, ny-4, 0, nz-1);

    boundaryCondition->addVelocityBoundary2N(water, *lattice, boundary::freeslip);
    boundaryCondition->addVelocityBoundary2P(top, *lattice, boundary::freeslip);
    boundaryCondition->addVelocityBoundary1N(rightWall, *lattice, boundary::freeslip);
    boundaryCondition->addVelocityBoundary1P(leftWall, *lattice, boundary::freeslip);

// //     boundaryCondition->setVelocityConditionOnBlockBoundaries(*lattice, top,       boundary::freeslip);
// //     boundaryCondition->setVelocityConditionOnBlockBoundaries(*lattice, rightWall, boundary::freeslip);
// //     boundaryCondition->setVelocityConditionOnBlockBoundaries(*lattice, leftWall,  boundary::freeslip);
    
    setBoundaryVelocity(*lattice, water,     Array<T,3>(uAveLB,(T)0.,(T)0.));
    setBoundaryVelocity(*lattice, top,       Array<T,3>(uAveLB,(T)0.,(T)0.));
    setBoundaryVelocity(*lattice, rightWall, Array<T,3>(uAveLB,(T)0.,(T)0.));
    setBoundaryVelocity(*lattice, leftWall,  Array<T,3>(uAveLB,(T)0.,(T)0.));
    
    integrateProcessingFunctional(new FluidPressureOutlet3D<T,DESCRIPTOR,0,+1>(), outlet1Domain, *lattice, 2);    
    setBoundaryVelocity(*lattice, outlet1Domain, Array<T,3>(uAveLB,(T)0.,(T)0.));
    
    
    iniLattice(*lattice, voxelizedDomain, uAveLB);
    if(iniVal) {
        Box3D toDomain(lattice->getBoundingBox());
        Box3D fromDomain(toDomain.shift(margin,margin,margin)); // During rescaling, the margin doubled in size,
                                                                //   an effect which is cancelled here through a shift.
        copy(*iniVal, fromDomain, *lattice, toDomain, modif::staticVariables);
    }


    
    // The ValueTracer is needed to check when a chosen quantity (in our case the average energy)
    //   has converged, so to conclude that steady state has been reached for the specific grid
    //   refinement level and stop the simulation.
    plint convergenceIter=20;
    util::ValueTracer<T> velocityTracer(0.05*convergenceIter, resolution, epsilon);
    global::timer("iteration").restart();
    plint i = util::roundToInt(currentTime/dt);
    lattice->resetTime(i);

    Box3D bbox(lattice->getBoundingBox());
    pcout << "Saving a " << lattice->getNx() << " by " << lattice->getNy()
          << " by " << lattice->getNz() << " lattice." << std::endl;
    global::timer("io").start();
    parallelIO::save(*lattice, "checkpoint", false);
    pcout << "Total time for i/o: " << global::timer("io").getTime() << std::endl;

    plint updateTimeBC = util::roundToInt(1/dt);
    plint outputTimeProfile = updateTimeBC;//util::roundToInt(1/dt);
    plint outputTimeVTK = updateTimeBC * 10;//util::roundToInt(10/dt);
    plint outputTimeFullVTK = updateTimeBC * 60;//util::roundToInt(60/dt);
    plint checkpointTime = updateTimeBC * 60;//util::roundToInt(60/dt);

    plint endTime = updateTimeBC * 600;// util::roundToInt(600/dt);

    // Collision and streaming iterations.
    pcout << "Starting "<< endTime <<" iterations" << std::endl;
    global::timer("global").start();

    for (plint i=0; i<endTime; ++i) {
		
                // if (fmod(i*dt,1)==0.) {
                if (fmod(i,updateTimeBC)==0.) { 
                    setBoundaryVelocity(*lattice, inletDomain, 
                        InletProfile(dx, dt));
                    pcout << "currentTime= " << currentTime << " updated inlet velocities" << endl;
                }
		        if (i%1000==0) {
                            pcout << "Write Image" << std::endl;
                            writeImages(*lattice,flagMatrix, lattice->getBoundingBox() , output+util::val2str(i), location, dx, dt );
                            pcout   << "; currentTime= " << currentTime
                                    // << "; av energy="
                                    // << setprecision(10) << getStoredAverageEnergy<T>(*lattice)
                                    << "; av rho="
                                    << setprecision(10) << getStoredAverageDensity<T>(*lattice) 
                                    << "; max. uLB= "
                                    << setprecision(10) << computeMax(*computeVelocityNorm(*lattice))
                                    << endl;
                        }
                // if (i%1000==0 && performOutput) {
                //             winddata(*lattice, location, i, dx, dt );
                            
                // }
                if (fmod(i,outputTimeProfile)==0.) {
                // if (i%5000000==0 && performOutput) {
                            windprofile(*lattice, location, i, dx, dt );
                }
                // if (i%5000000==0 && performOutput) {
                if (fmod(i,outputTimeVTK)==0.) {
                            pcout << "Write VTK" << std::endl;
                            // writeVTK(*lattice,flagMatrix, lattice->getBoundingBox() , output+util::val2str(i), location, dx, dt );
                            writeVTK(*lattice,flagMatrix, output+util::val2str(i), location, dx, dt );
                }
                if (fmod(i,outputTimeFullVTK)==0.) {
                // if (i%5000000==0 && performOutput) {
                            pcout << "Write FullVTK" << std::endl;
                            writeFullVTK(*lattice,flagMatrix, lattice->getBoundingBox() , output+"_Full_"+util::val2str(i), location, dx, dt );
                }
                if (fmod(i,checkpointTime)==0. && i>1) { 
                           pcout << "Saving the state of the simulation ..." << endl;
            //saveRawMultiBlock(lattice, "checkpoint.dat");
                    saveBinaryBlock(*lattice, "checkpoint.dat"); 
                }
        lattice->collideAndStream();

        ++i;
        currentTime = i*dt;
    }
    pcout << "End of simulation" << std::endl;
    pcout << "Total time of execution: " << global::timer("global").getTime() << std::endl;

    // Box3D measureBox(inletDomain);

    return lattice;
}
	
// Read the user input XML file provided at the command-line.
void readParameters(XMLreader const& document)
{
    std::string meshFileName;
    std::vector<std::string> openingType;
    document["geometry"]["mesh"].read(meshFileName);
    document["geometry"]["inletDiameter"].read(userDefinedInletDiameter);
    document["geometry"]["averageInletVelocity"].read(averageInletVelocity);
    document["geometry"]["openings"]["sortDirection"].read(openingSortDirection);
    document["geometry"]["openings"]["type"].read(openingType);

    document["fluid"]["kinematicViscosity"].read(kinematicViscosity);
    document["fluid"]["density"].read(fluidDensity);
    document["fluid"]["volume"].read(volume);

    document["numerics"]["referenceDirection"].read(referenceDirection);
    document["numerics"]["referenceResolution"].read(referenceResolution);
    document["numerics"]["nuLB"].read(nuLB);

    document["simulation"]["simTime"].read(simTime);
    document["simulation"]["maxLevel"].read(maxLevel);
    document["simulation"]["epsilon"].read(epsilon);

    document["simulation"]["performOutput"].read(performOutput);
    document["simulation"]["doImages"].read(doImages);
    document["simulation"]["useAllDirections"].read(useAllDirections);
    document["simulation"]["useRegularizedWall"].read(useRegularizedWall);
    document["simulation"]["useIncompressible"].read(useIncompressible);
    document["simulation"]["poiseuilleInlet"].read(poiseuilleInlet);
    document["simulation"]["convectiveScaling"].read(convectiveScaling);

    // At this part, the surface geometry of the aneurysm (as given by the user in
    //   the form of an ASCII or binary STL file) is read into a data structure
    //   comprised by a set of triangles. The DBL constant means that double
    //   precision accuracy will be used (generally the recommended choice).
    triangleSet = new TriangleSet<T>(meshFileName, DBL);
    pcout << "Reynolds number, based on provided inlet diameter: "
          << averageInletVelocity*userDefinedInletDiameter/kinematicViscosity
          << std::endl;
    plbIOError(openingSortDirection<0 || openingSortDirection>2,
               "Sort-direction of opening must be 0 (x), 1 (y), or 2 (z).");
    // The surface geometry, as provided by the STL file, must contain openings,
    //   namely inlets and outlets. On these openings, appropriate boundary conditions
    //   will be imposed by palabos. Which opening is inlet and which is outlet, is
    //   identified by the user in the input XML file.
    openings.resize(openingType.size());
    for (pluint i=0; i<openingType.size(); ++i) {
        std::string next_opening = util::tolower(openingType[i]);
        if (next_opening=="inlet") {
            openings[i].inlet = true;
        }
        else if (next_opening=="outlet") {
            openings[i].inlet = false;
        }
        else {
            plbIOError("Unknown opening type.");
        }
    }
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::IOpolicy().activateParallelIO(true);

    

    string paramXmlFileName;
    string outFolder;
    try {
        global::argv(1).read(outFolder);
        global::argv(2).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " output_older parameter-input-file.xml" << std::endl;
        return -1;
    }
    
    std::string outputDirectory = outFolder+"/";
    global::directories().setOutputDir(outputDirectory.c_str());
    
    
    
    // Read the parameter XML input file. (Lots of comments are included there too).
    try {
        XMLreader document(paramXmlFileName);
        readParameters(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Error in input file " << paramXmlFileName
              << ": " << exception.what() << std::endl;
        return -1;
    }

    plint iniLevel=0;
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > iniConditionLattice(0);
    // This code incorporates the concept of smooth grid refinement until convergence is
    //   achieved. The word ``smooth'' indicates that as the refinement level increases
    //   by one, the whole grid doubles in each direction. When the grid is refined, both
    //   dx and dt have to change. Whether dt is changed as dx^2 (diffusive behavior)
    //   or as dx (convective behavior), is controlled by the input variable
    //   ``convectiveScaling'' (the recommended choice is not to use convective scaling).
    try {
        for (plint level=iniLevel; level<=maxLevel; ++level) {
            pcout << std::endl << "Running new simulation at level " << level << std::endl;
            std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > convergedLattice (
                    run(level, iniConditionLattice.get()) );
            if (level != maxLevel) {
                plint dxScale = -1;
                plint dtScale = -2;
                if (convectiveScaling) {
                    dtScale = -1;
                }
                // The converged simulation of the previous grid level is used as the initial condition
                //   for the simulation at the next grid level (after appropriate interpolation has
                //   taken place).
                iniConditionLattice = std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > (
                        refine(*convergedLattice, dxScale, dtScale, new BGKdynamics<T,DESCRIPTOR>(1.)) );
            }
        }
    }
    catch(PlbException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}

