#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sstream>

using namespace std;

int main(int argc, char ** argv)
{

  if(argc < 6)
  {
    cout << "Usage: exec topfile numresidue numatoms numframes cutoff [outfile]" << endl;
    return 1;
  }

  int rsize = atoi(argv[2]);
  int asize = atoi(argv[3]);
  int fsize = atoi(argv[4]);
  double ctoff = atof(argv[5]);
  char pdbfile[70];
  char ofile[70];

  if(argc > 6) {
     strcpy(ofile, argv[6]);
  } else {
     strcpy(ofile, "contact.trajectory");
  }

  strcpy(pdbfile, argv[1]);

  ifstream in(pdbfile);

  char dline[300];
  int rid[asize];
  int acount = 0;
  int rcount = 0;
  while (true)
  {
    in.getline(dline, 300);
    stringstream instr;
    instr.str(dline);

    string trash;

    instr >> trash;
    if (trash == "ATOM")
    {
      for (int k = 0; k < 3; k++)
      {
        instr >> trash;
      }

      instr >> rid[acount];
      if(acount==0 || rid[acount] != rid[acount-1]) {
         rcount = rcount + 1; 
      }
      acount = acount + 1;
      //for (int k = 0; k < 3; k++)
      //{
      // instr >> trash;
      //}
    }

    if (trash == "END" || instr.eof())
    {
      break;
    }
  }
  in.close();

//  cout << rcount << acount << endl; 

  if ( asize != acount )
  {
    cout << " Atom number mismatch  " << endl;
    return 1;
  }
  if ( rsize != rcount )
  {
    cout << " Residue number mismatch  " << endl;
    return 1;
  }


  //int csize = rsize*(rsize + 1)/2;
  int mcontacts[rsize][rsize];

  ifstream inm("meaningful.contacts");
  if (!inm.is_open())
  {
    cout << "ERROR: Can't find meaningful.contacts file! " << endl;
    return 1;
  }

  for(int i = 0; i < rsize; i++)
  {
    for(int j = 0; j < rsize; j++)
    {
      mcontacts[i][j] = 0;
    }
  }

  while(true) //for(int i=0; i < csize; i++)
  {
    int iindex = 0;
    int jindex = 0;
    inm >> iindex;

    if(inm.eof())
    {
      break;
    }

    inm >> jindex;

    mcontacts[iindex - 1][jindex - 1] = 1;
  }

  inm.close();

//  ofstream out("contact.trajectory");
  ofstream out(ofile);

  double x[asize];
  double y[asize];
  double z[asize];

  for(int i0 = 0; i0 < fsize; i0++) 
  {
    for(int j0 = 0; j0 < asize; j0++) 
    {
      cin >> x[j0] >> y[j0] >> z[j0];
//      x[j0] = round(x[j0]*1000) / 1000;
//      y[j0] = round(y[j0]*1000) / 1000;
//      z[j0] = round(z[j0]*1000) / 1000;
    }
    int pdbcontact[rsize][rsize];
    for(int j = 0; j < rsize; j++)
    {
      for(int k = 0; k < rsize; k++)
      {
          pdbcontact[j][k] = 0;
      }
    }

    for(int j = 0; j < asize; j++)   //check for contact
    {
      for(int k = j; k < asize; k++) // no redundant
      {
        if(mcontacts[rid[j] - 1][rid[k] - 1] == 1)
        {
          if (((x[j]-x[k])*(x[j]-x[k]) + 
               (y[j]-y[k])*(y[j]-y[k]) + 
               (z[j]-z[k])*(z[j]-z[k])) < ctoff*ctoff) 
          {
            pdbcontact[rid[j] - 1][rid[k] - 1] = 1;
          }
        }
      }
    }

    for(int j = 0; j < rsize; j++)   //check for contact
    {
      for(int k = j; k < rsize; k++) // no redundant
      {
        if (mcontacts[j][k] == 1)
        {
          out << pdbcontact[j][k] << " "; 
        }
      }
    }
    out << endl;
  }

  out.close();

  return 0;
}
