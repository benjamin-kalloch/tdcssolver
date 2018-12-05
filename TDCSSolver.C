#include "fvCFD.H"
#include "simpleControl.H"
// this include is neccessary so 'transform()' will operate on the GeometricFields (= volVectorField) not the super-class Field-type
// otherwise the return value would be a Field as well and thus cannot be assigned to J again (which in turn is a GeometricField...)
#include "transformGeometricField.H"

#include <algorithm>
#include <vector>
#include <memory>
#include <numeric>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main( int argc, char *argv[] )
{
    #include "setRootCase.H"

    #include "createTime.H"     // Initializes the global 'runTime'-object controlling all timestep related things
    #include "createMesh.H"     // Initializes the global 'mesh'-object controlling all mesh related things
    #include "createFields.H"   // Custom header file to intialize the necessarry fields

    // Simple control is a subclass of solutionControl implementing the SIMPLE algorithm
    // However, we only use the capabilities of solutionControl which must be instanciated
    // by providing the name of a specific solution control algorithm;
    simpleControl simple(mesh);


    Info<< "Calculating distribution of electrical potential." << nl << endl;
    while( simple.loop( runTime ) )
    {
        Info<< "Iteration = " << runTime.timeName() << endl;

        while( simple.correctNonOrthogonal() )
        {
            fvScalarMatrix ElPotEqn
            (
                isScalarField ? 
                    fvm::laplacian( *(std::static_pointer_cast< volScalarField >(sigma)), ElPot) 
                    : 
                    fvm::laplacian( *(std::static_pointer_cast< volTensorField >(sigma)), ElPot)
            );
            // The previous state of the field must  be saved prior to solving for successful relaxation
            ElPot.storePrevIter();
            ElPotEqn.solve();

            // use under-relaxation to stabilize the result 1=mild relexation, 0=heavy relaxation,
            // The factor can be hardcoded (as parameter to relax()) or as in this case read from 
            // relaxationFactors-dictonary in fvSolution if not corresponding dictionary is provided,
            //  the relaxation will be skipped
            ElPot.relax();
        }
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
            << " ClockTime = "    << runTime.elapsedClockTime() << " s"
            << nl << endl;
#ifdef DEBUG
        E = - fvc::grad( ElPot );
        ElPot.write();
        E.write();
#endif
    }
 
    Info<< "Calculating gradients... ";
    E = - ( fvc::grad( ElPot) );             // gradient field of the smoothed ElPot
    if( isScalarField )
        J = *(std::static_pointer_cast< volScalarField >(sigma)) * E;
    else
        J = transform( *(std::static_pointer_cast< volTensorField >(sigma)), E );
    Info<< "Done." << nl << endl;

    Info<< "Determining scaling factor." << endl;

    // Volume integral over entire domain
    Foam::dimensioned< Foam::Vector<scalar> > volIntegral = fvc::domainIntegrate( J );
    Info << "Volume integral of J over entire domain: " << volIntegral.value() << endl;

    // Integral over the electrode-surfaces
    scalar jIntegrated = numElectrodes > 0 ? 0. : targetCurrent;    // scale by 1 if no surface to integrate over are provided
   
    // Interpolate values of the control volumes onto their faces
    surfaceVectorField sField = fvc::interpolate(J);

    // Integral over the faces representing the contact surface of the electrodes
    scalar perElectrodeIntegrated = 0;
    for(size_t e = 0; e < numElectrodes; e++)
    {
        label faceZoneID = mesh.faceZones().findZoneID( patchNames[e] );
        const labelList &faces = mesh.faceZones()[faceZoneID];    

        std::vector<scalar> faceValues( faces.size() );
        forAll( faces, f )
        {
            label faceI = faces[ f ];        // & = dot product
            faceValues[f] = fabs( sField[faceI] & mesh.Sf()[faceI] );
        }
        
        std::cout << "Begin: " << faceValues[0] << std::endl;

        perElectrodeIntegrated = std::accumulate(faceValues.begin(), faceValues.end(), scalar(0.0) );

        Info<< "Integral over faceZone (ID: " << faceZoneID <<" ) '"<< patchNames[e]
            << "' =" << perElectrodeIntegrated << endl;
        
        jIntegrated += perElectrodeIntegrated;
    }

    // Determine the factor between the integrated current density at 
    // the electrode surface and the target current density that should
    // be applied through the electrodes.
    // We use this value to scale the entire field accordingly.  
    scalar factor = targetCurrent / jIntegrated;

    Info<< nl << "Scaling current density field by factor "<< factor << nl << endl;

    J *= factor;
    E *= factor;

    magE = mag(E);
    magJ = mag(J);

    // Print runtime statistics
    Info<< "**** Simulation finshed! ****" << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime()   << " s"
        << " ClockTime = "    << runTime.elapsedClockTime() << " s"
        << nl << endl;
    
    Info<< "Writing results." << nl << endl;
    runTime++;
    J.write();
    magJ.write();
    E.write();
    magE.write();
    ElPot.write();

    if( isScalarField )
        std::static_pointer_cast< volScalarField >(sigma)->write();
    else
        std::static_pointer_cast< volTensorField >(sigma)->write();
   

    return(0);
}
// ************************************************************************* //
