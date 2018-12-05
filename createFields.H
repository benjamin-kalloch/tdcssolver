    // *** Handling the solver properties ***
    IOdictionary solverProperties
    (
        // further information on 'IOobject':
        // https://openfoamwiki.net/index.php/OpenFOAM_guide/Input_and_Output_operations_using_dictionaries_and_the_IOobject_class
        IOobject
        (
            "solverProperties",     // name of the file containing the dictionary
            runTime.system(),       // the location of the dictionary file; here: the constant-directory
            mesh,                   // the objectRegistry is represented by the mesh
            IOobject::MUST_READ,    // the object must be read Istream on construction - MUST_READ_IF_MODIFIED is also possible
            IOobject::NO_WRITE      // the object is not written on destruction, but can be written explicitly
        )
    );
    
    scalar targetCurrent  = readScalar( solverProperties.lookup("targetCurrent") );
    size_t numElectrodes  = static_cast<size_t>(readScalar( solverProperties.lookup("numElectrodes") ));
    Info<<"Using solver properties:\n\t> input urrent\t: " << targetCurrent << " ampere\n"
        <<"\t> #electrodes\t: " << numElectrodes << "\n"
        <<"\t> patch names\t: " << endl;

    word patchNames[numElectrodes];
    ITstream &patches =  solverProperties.lookup("electrodeContactSurfaceNames");          
    token temp;
    for(size_t w = 0; w < numElectrodes; w++)
    {
        
        patches.read(temp);
        patchNames[w] = temp.wordToken();
        Info<<"\t\t- " << patchNames[w] << endl;
    }

    // *** Handling field sigma ***
    std::shared_ptr<void> sigma;
    IOobject sigma_IO
    (
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );
    
    // determine type of sigma
    bool isScalarField = true;
    
    Info<< "Reading field sigma as '";
    std::ifstream sigma_tmp( sigma_IO.objectPath() );
    std::string line;
    while( std::getline( sigma_tmp, line) )
    {
        if( line.find( "volScalarField" ) != std::string::npos )
        {
            Info << "volScalarField'." << endl;
            break;
        }
        else if ( line.find( "volTensorField" ) != std::string::npos)
        {
            Info << "volTensorField'." << endl;
            isScalarField = false;
            break;
        }
    }

    if( isScalarField)
    {
        sigma = std::static_pointer_cast<void>(std::make_shared<volScalarField>(
            sigma_IO,
            mesh
        ));
    }
    else
    {
        sigma = std::static_pointer_cast<void>(std::make_shared<volTensorField>(
            sigma_IO,
            mesh
        ));
    }

    // *** Handling field ElPot ***
    Info<< "Reading field ElPot\n" << endl;
    volScalarField ElPot
    (
        IOobject
        (
            "ElPot",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    
    // *** Handling field E ***
    dimensionedVector E_units( "zero", dimensionSet(1, 1, -3, 0, 0, -1, 0),vector::zero );
    dimensionedScalar magE_units( "zero", dimensionSet(1, 1, -3, 0, 0, -1, 0), 0 );
    
    volVectorField E
    (
        IOobject
        (
            "E",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        E_units
    );

    volScalarField magE
    (
        IOobject
        (
            "magE",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        magE_units
    );

    // *** Handling field J ***
    dimensionedVector J_units( "zero", dimensionSet(0, -2, 0, 0, 0, 1, 0),vector::zero );
    dimensionedScalar magJ_units( "zero", dimensionSet(0, -2, 0, 0, 0, 1, 0), 0 );
    
    volVectorField J
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        // When representing sigma as a tensor, we must use the same dimension as the electrical field strength.
        // Because in this case we transform the E field by sigma and do not multiply which would change the units of measurement.
        (isScalarField) ? J_units : E_units
    );

    volScalarField magJ
    (
        IOobject
        (
            "magJ",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        (isScalarField) ? magJ_units : magE_units
    );
