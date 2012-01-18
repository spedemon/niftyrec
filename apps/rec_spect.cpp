
#include <string>
#include <iostream>
#include <algorithm>
#include <tclap/CmdLine.h>
#include "niftyrec_version.h"
#include "_et_array_interface.h"

#ifdef _SUPPORT_NRRD
#include "teem/nrrd.h"
#endif 

#define FILETYPE_NIFTI 0
#define FILETYPE_NRRD  1

#define XML_FILE "./rec_spect.xml"

char xml_string[] = ""
"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
"<executable>\n"
"  <category>Reconstruction</category>\n"
"  <title>REC SPECT</title>\n"
"  <description>\n"
"  Reconstruct SPECT images from sinogram data\n"
"  </description>\n"
"  <version>1.4</version>\n"
"  <documentation-url>http://sourceforge.net/niftyrec</documentation-url>\n"
"  <license>http://sourceforge.net/niftyrec</license>\n"
"  <contributor>UCL</contributor>\n"
"  <parameters>\n"
"    <label>Workspace IO</label>\n"
"    <description>Input images</description>\n"
"    <image>\n"
"      <name>Sinogram</name>\n"
"      <label>Sinogram</label>\n"
"      <channel>input</channel>\n"
"      <longflag>sinogram</longflag>\n"
"      <description>Sinogram data</description>\n"
"    </image>\n"
"    <image>\n"
"      <name>OutputActivity</name>\n"
"      <label>Output Activity</label>\n"
"      <channel>output</channel>\n"
"      <longflag>output</longflag>\n"
"      <description>Output activity</description>\n"
"    </image>\n"
"  </parameters>\n"
"  <parameters>\n"
"    <label>SPECT Acquisition Parameters</label>\n"
"    <description>\n"
"    Parameters of the SPECT acquisition\n"
"    </description>\n"
"    <double>\n"
"      <name>FirstCamera</name>\n"
"      <longflag>firstcamera</longflag>\n"
"      <description>\n"
"      First camera position in degrees\n"
"      </description>\n"
"      <label>First Camera Position</label>\n"
"      <default>0.0</default>\n"
"    </double>\n"
"    <double>\n"
"      <name>LastCamera</name>\n"
"      <longflag>lastcamera</longflag>\n"
"      <description>\n"
"      Last camera position in degrees\n"
"      </description>\n"
"      <label>Last Camera Position</label>\n"
"      <default>180.0</default>\n"
"    </double>\n"
"  </parameters>\n"
"  <parameters>\n"
"    <label>Reconstruction Parameters</label>\n"
"    <description>\n"
"    Parameters for the reconstruction\n"
"    </description>\n"
"    <integer>\n"
"      <name>iterations</name>\n"
"      <longflag>iterations</longflag>\n"
"      <description>\n"
"      An integer without constraints\n"
"      </description>\n"
"      <label>Iterations</label>\n"
"      <default>20</default>\n"
"    </integer>\n"
"  </parameters>\n"
"</executable>\n";




int stream_xml()
{
        //if the XML file is in the path, read it and stream it to stdout, otherwise use default from compile time
        std::string line;
        std::ifstream xml_file (XML_FILE);
        if (xml_file.is_open())
        {
            while ( xml_file.good() )
            {
                std::getline (xml_file,line);
                std::cout << line << std::endl;
            }
        xml_file.close();
        }
        else 
        {
        std::cout << xml_string;
        }
        return 0;   
}

nifti_image *nrrd_to_nifti(Nrrd *nrrdImage)
{
    nifti_image *niftiImage = NULL;
    /* say something about the array */
    std::cout << "== Nrrd -> Nifti conversion ==" << std::endl;
    printf("-- Converting %d-dimensional NRRD of type %d (%s)\n", nrrdImage->dim, nrrdImage->type, airEnumStr(nrrdType, nrrdImage->type));
    for (int n_axis=0; n_axis<nrrdImage->dim; n_axis++)
        printf("-- Size axis %d: %d\n",n_axis,nrrdImage->axis[n_axis].size);
    printf("-- The array contains %d elements, each %d bytes in size\n", (int)nrrdElementNumber(nrrdImage), (int)nrrdElementSize(nrrdImage));
    /* convert it: size and dimensions */
    int dim[8];
    dim[0] = nrrdImage->dim;
    for (int n_axis=0; n_axis<nrrdImage->dim; n_axis++)
        dim[n_axis+1]=nrrdImage->axis[n_axis].size;
    for (int n_axis=nrrdImage->dim+1; n_axis<8; n_axis++)
        dim[n_axis]=1;
    /* convert it: spatial transformations */
    //.. set axis.measurementFrame, axis.spaceOrigin, axis.spaceUnits
    if (nrrdImage->type == nrrdTypeFloat) 
    {
        //allocate new Nifti image struct
        niftiImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        //copy data (rather than passing the pointer, to simplify free - use nrrdNuke() and nifti_image_free() )
        niftiImage->data = (void *) calloc((int)nrrdElementNumber(nrrdImage),sizeof(float));
        memcpy((void*)niftiImage->data, (void*)nrrdImage->data,(int)nrrdElementNumber(nrrdImage)*sizeof(float));
        std::cout << "-- Conversion done" << std::endl;
    }
    else
    {
        std::cout << "-- Conversion from Nrrd to Nifti not supported for Nrrd format %s." << airEnumStr(nrrdType, nrrdImage->type) << std::endl;
    }
    return niftiImage;
}


Nrrd *nifti_to_nrrd(nifti_image *niftiImage)
{
    std::cout << "== Nifti -> Nrrd conversion ==" << std::endl;
    Nrrd *nrrdImage = nrrdNew(); 
    unsigned int dim = niftiImage->dim[0];
    size_t size[NRRD_DIM_MAX];
    for (int n_axis=0; n_axis<dim; n_axis++)
        size[n_axis] = niftiImage->dim[n_axis+1];
    if(niftiImage->datatype == NIFTI_TYPE_FLOAT32)
    {
        if (nrrdAlloc_nva(nrrdImage, nrrdTypeFloat, dim, size)) 
        {
            std::cout << "-- Conversion failed, cannot allocate the NRRD image. " << std::endl;
            return NULL;
        } 
        /* copy data - rather than passing the pointer, to facilitate a bit memory handling */
        memcpy((void*)nrrdImage->data,(void*)niftiImage->data,(int)nrrdElementNumber(nrrdImage)*sizeof(float));
    }
    else
    {
        std::cout << "-- Conversion from Nifti to Nrrd not supported for Nifti data type %d." << niftiImage->datatype << std::endl;
        return NULL;
    }
    printf("-- Created %d-dimensional NRRD of type %d (%s)\n", nrrdImage->dim, nrrdImage->type, airEnumStr(nrrdType, nrrdImage->type));
    for (int n_axis=0; n_axis<nrrdImage->dim; n_axis++)
        printf("-- Size axis %d: %d\n",n_axis,nrrdImage->axis[n_axis].size);
    std::cout << "-- Conversion done" << std::endl;
    return nrrdImage;
}





int main(int argc, char** argv)
{
    if (argc==2)
    {
        if(strcmp(argv[1], "--xml") == 0)
        {
            stream_xml();
            return 0;
        }
    }       
        
    try 
    {  
    //Parse command line arguments
        TCLAP::CmdLine cmd("Maximum Likelihood Expectation Maximisation (MLEM) reconstruction. ", ' ', VERSION);

        //Required input sinogram data file
        TCLAP::ValueArg<std::string> sinogramArg("s","sinogram","Sinogram data file to reconstruct",true,"homer","string",cmd);
        //Required output file name
        TCLAP::ValueArg<std::string> outputArg("o","output","File name for the output reconstructed image",true,"homer","string",cmd);
        //Optional first camera angle
        TCLAP::ValueArg<float> firstcameraArg("f","firstcamera","Position of the first camera in degrees. Default is 0.0 deg",false,0.0f,"float",cmd);
        //Optional last camera angle
        TCLAP::ValueArg<float> lastcameraArg("l","lastcamera","Position of the last camera in degrees. Default is 180.0 deg",false,180.0f,"float",cmd);
        //Optional number of iterations
        TCLAP::ValueArg<int> iterationsArg("i","iterations","MLEM iterations",false,10,"int",cmd);
        //Optional switch verbose
        TCLAP::SwitchArg verboseSwitch("v","verbose","Print verbose", cmd, false);
        //Optional switch use gpu
        TCLAP::SwitchArg gpuSwitch("g","gpu","Use Graphics Processing Unit", cmd, false);


        cmd.parse( argc, argv );

        std::string sinogram_filename = sinogramArg.getValue().c_str();
        std::string output_filename = outputArg.getValue().c_str();
        float firstcamera = firstcameraArg.getValue();
        float lastcamera = lastcameraArg.getValue();
        int iterations = iterationsArg.getValue();
        bool verbose = verboseSwitch.getValue();


        if (verbose)
        {
            std::cout << "Sinogram file name: " << sinogram_filename << std::endl;
            std::cout << "Output file name:   " << output_filename << std::endl;
            std::cout << "First camera:       " << firstcamera << std::endl;
            std::cout << "Last camera:        " << lastcamera << std::endl;
            std::cout << "Iterations:         " << iterations << std::endl;
        }

        //Load input sinogram
        int filetype_sinogram;
        nifti_image *sinogramImage = NULL;

        if (sinogram_filename.substr(sinogram_filename.find_last_of(".") + 1) == "nii" || sinogram_filename.substr(sinogram_filename.find_last_of(".") + 1) == "nii.gz")
        {
            std::cout << "Nifti sinogram file..." << std::endl;
            filetype_sinogram = FILETYPE_NIFTI;
        }
        else if (sinogram_filename.substr(sinogram_filename.find_last_of(".") + 1) == "nrrd")
        {
            std::cout << "NRRD sinogram file..." << std::endl;
            filetype_sinogram = FILETYPE_NRRD;
        }
        else
        {
            std::cout << "Unknown file format" << std::endl;
            return 1;
        }

        if (filetype_sinogram == FILETYPE_NIFTI)
        {
            sinogramImage = nifti_image_read(sinogram_filename.c_str(),true);
            if (sinogramImage == NULL)
            {
                std::cout << "Couldn't load the sinogram " << sinogram_filename << std::endl;
                return 1;                
            }
        }
        else if (filetype_sinogram == FILETYPE_NRRD)
        {
#ifdef _SUPPORT_NRRD
            /* create a nrrd; at this point this is just an empty container */
            Nrrd *sinogram_nrrdImage;
            sinogram_nrrdImage = nrrdNew();
            char *err;

            /* read in the nrrd from file */
            if (nrrdLoad(sinogram_nrrdImage, sinogram_filename.c_str(), NULL))
            {
                err = biffGetDone(NRRD);
                fprintf(stderr, "Problem reading \"%s\":\n%s", sinogram_filename.c_str(), err);
                free(err);
                return 1;
            }

            /* convert NRRD to Nifti */

            sinogramImage = nrrd_to_nifti(sinogram_nrrdImage);
            if (sinogramImage == NULL)
            {
                std::cout << "Convertion from Nrrd to Nifti failed." << std::endl;
                return 1;
            }

            /* blow away both the Nrrd struct *and* the memory at nin->data. (nrrdNix() frees the struct but not the data, nrrdEmpty() frees the data but not the struct) */
            nrrdNuke(sinogram_nrrdImage);
#else 
            std::cout << "NRRD file format not supported. Please enable the NRRD support flag at compile time." << std::endl;
            return 1;
#endif 
        }
        else
        {
            std::cout << "Unknown file format" << std::endl;
            return 1;            
        }


        //Reconstruct
        float *sinogram_data = (float*)sinogramImage->data;
        int size_x = sinogramImage->nx;
        int size_y = sinogramImage->ny;
        int n_cameras = sinogramImage->nz;
        int use_psf = 0;
        int use_ddpsf = 1;
        int psf_size_x; 
        int psf_size_y;
        float *psf_data;
        int use_attenuation = 0; 
        float *attenuation_data;
        int use_gpu = (int)verboseSwitch.getValue();;

        printf("== MLEM Reconstruction == \n");
        printf("-- size_x: %d   \n-- size_y: %d   \n-- n_cameras: %d   \n-- use_psf: %d   \n", size_x, size_y, n_cameras, use_psf);
        printf("-- use_ddpsf: %d   \n-- use_attenuation: %d   \n-- use_gpu: %d    \n", use_ddpsf, use_attenuation, use_gpu);

        float *activity_data = (float*) calloc(size_x*size_x*size_y,sizeof(float));
        if (et_mlem_spect(sinogram_data, size_x ,size_y ,n_cameras, firstcamera, lastcamera, iterations, use_psf, use_ddpsf, psf_size_x, psf_size_y, psf_data, use_attenuation, attenuation_data, activity_data, use_gpu))
        {
            std::cout << "MLEM reconstruction failed. " << std::endl;
            return 1;
        }

        //Save result
        int filetype_output;

        if (output_filename.substr(output_filename.find_last_of(".") + 1) == "nii" || output_filename.substr(output_filename.find_last_of(".") + 1) == "nii.gz")
        {
            std::cout << "Nifti output file ..." << std::endl;
            filetype_output = FILETYPE_NIFTI;
        }
        else if (output_filename.substr(output_filename.find_last_of(".") + 1) == "nrrd")
        {
            std::cout << "NRRD output file ..." << std::endl;
            filetype_output = FILETYPE_NRRD;
        }
        else
        {
            std::cout << "Unknown file format" << std::endl;
            return 1;
        }

        int dim[8];
        dim[0]=3;
        dim[1]=size_x;
        dim[2]=size_x;
        dim[3]=size_y;
        dim[4]=1; dim[5]=1; dim[6]=1; dim[7]=1;
        nifti_image *activityImage;
        activityImage = nifti_make_new_nim(dim, NIFTI_TYPE_FLOAT32, false);
        activityImage->data = static_cast<void *>(activity_data);        
        nifti_set_filenames(activityImage, output_filename.c_str(), 0, 0);

        if (filetype_output == FILETYPE_NIFTI)
        {
            nifti_image_write(activityImage);
        }
        else if (filetype_output == FILETYPE_NRRD)
        {
#ifdef _SUPPORT_NRRD
            /* write out the nrrd to a different file */
            Nrrd *activity_nrrdImage = nifti_to_nrrd(activityImage);
            char *err;
//            nrrdLoad(activity_nrrdImage, sinogram_filename.c_str(), NULL);//
            if (nrrdSave(output_filename.c_str(), activity_nrrdImage, NULL))
            {
                err = biffGetDone(NRRD);
                fprintf(stderr, "Trouble writing \"%s\":\n%s", output_filename.c_str(), err);
                free(err);
                return 1;
            }
            nrrdNuke(activity_nrrdImage);
#else
            std::cout << "NRRD file format not supported. Please enable the NRRD support flag at compile time." << std::endl;
            return 1;
#endif
        }
        //Clean up
        nifti_image_free(activityImage);

    }
    catch (TCLAP::ArgException &e)  // catch any exceptions
    {   std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; } 
    return 0;
}


