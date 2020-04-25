
/* Finding the minimal basis of a single linear Diophantine equation
   of the form  sum(ai xi) = sum(bj yj)
   by the APT2 algorithm
*/
/*---------------------------------------------------------------------------*\

   Copyright (c) 1991--2000 Ana Paula Tomas, Miguel Filgueiras /
  DCC-FC & LIACC, Universidade do Porto

#    This program is free software; you can redistribute it and/or modify
#      it under the terms of the GNU General Public License as published by
#      the Free Software Foundation; either version 2 of the License, or
#      (at your option) any later version.
#
#      This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#      GNU General Public License for more details.
#
#      You should have received a copy of the GNU General Public License
#      along with this program.

\*---------------------------------------------------------------------------*/

/* The APT2 algorithm is
   the super-generalized APT algorithm with slope information,
   slopes are computed in axbjyj, also when v<0
   - with Lambert bound everywhere and using bounds for y in aptfilter

   NEW VERSION: solving "k.dy < y mod ymax" to get the appropriate slope

   - with array of comparable solutions in deadfilter/deadsfilter

   - using the "good order" for the coefficients
*/

#define VERSION "Slopes-V7i"

#include <stdio.h>

/* Comment/uncomment next line to not count/count candidates */
/* #define COUNTCAND 1 */

#define MAXDIM 32
#define MAXSOL 50000
#define MAXSLOP 100000

#define ALL (MAXSOL+2)
#define NONE (-1)

#define ISDIGIT(C) (C>='0' && C<='9')
#define MODASS(X,Y) if ((X=X%Y) < 0)  X+=Y

typedef unsigned long BITS;

typedef struct el {
  int cmax;  BITS cmask;  struct el *nxt;
} EL, *LIST;

int Coefs[MAXDIM], N, Ni, Nj, Dead,
  CVars, ConfVars[MAXDIM], ConfVals[MAXDIM],
  Sols[MAXSOL][MAXDIM], NSols, L2Sol, CLSol, CCSol,
  Max[MAXDIM], GCD[MAXDIM][MAXDIM],
  GCDa2[MAXDIM][MAXDIM], MultA2[MAXDIM][MAXDIM], Mult[MAXDIM][MAXDIM],
  Dz[MAXDIM][MAXDIM], Dy[MAXDIM][MAXDIM], YSlop[MAXSLOP], ZSlop[MAXSLOP],
  VSlop[MAXSLOP], VLgDz, SmDy, 
  E0, ENj, ERhsL, EDead, Ea, Eb, Ec, Eix, EMax[MAXDIM], EL2Sol,
  XMaxi, XMaxj, XGCDai, XGCDa2ij, XDzij, XDyij, XMultAi, XMultA2ij,
  XGCDAj, XMultAj,
  CMaxi, CMaxj, CMaxs[MAXDIM], LhsL, LhsSum,
  LhsSV, RhsSV, LhsSup, RhsSup, ESV,
  CompSols[MAXSOL], LastCmpSol, NextCmpSol,
  InOrder[MAXDIM], OutOrder;

BITS NZSol[MAXSOL], DeadConfs[MAXSOL], Mask, Bit[MAXDIM], EDeadCs[MAXSOL];

LIST LMaxs[MAXDIM], ELMaxs[MAXDIM];

char Version[] = VERSION;

#ifdef COUNTCAND
double Rej0, RejF;
#endif

LIST newel(mx,msk,n)
int mx;  BITS msk;  LIST n;
{ LIST l = (LIST) malloc(sizeof(EL));

  l->cmax = mx;  l->cmask = msk;  l->nxt = n;
  return(l);
}

int inint()
{ int k=0, p;
  char c;

  do c=getchar();
  while ((p=(c!='-')) && ! ISDIGIT(c));
  if (!p) {
    c=getchar();
    if (! ISDIGIT(c)) {
      printf("Bad char after - sign: %c\n",c);
      exit(0);
    }
  }
  do {
    k=c-'0'+k*10;  c=getchar();
  } while (ISDIGIT(c));
  if (p) return(k);
  else return(-k);
}

void inputdata()
{ int *pp = Coefs;

  Ni = Nj = 0;
  while ((*pp++ = inint()))  Ni++;
  pp--;
  while ((*pp++ = inint()))  Nj++;
}

void printconf(i,j)
int i, j;
{ printf(">>configuration: y%d, y%d\tMask=%o\n\t",i+1,j+1,Mask);
  for (i=0; i<CVars; i++)
    printf("y%d=%d ",ConfVars[i]+1,ConfVals[i]);
  putchar('\n');
}

int actual(i)
int i;
{ int k;

  if (OutOrder) {
    for (k=0; InOrder[k]!=i; k++);
    return(k);
  }
  return(i);
}
    
void printsol(k)
int k;
{ int i, ii;

  printf("( ");
  for (i=0; i<N; i++) {
    ii = actual(i);
    if (Bit[ii]&NZSol[k]) printf("%d ",Sols[k][ii]);
    else printf("0 ");
  }
  printf(")\n");
}

void ordercoeffs()
{ int i = LhsSup, j, k, m, c;

  if (i < RhsSup)  i = RhsSup;
  if (i>100 || N>6 && i>20) {  /* force "good order" */
    if (!OutOrder) {   /* this is only to avoid problems when repeating */
      OutOrder = 1;
      for (i=0; i<N; i++)  InOrder[i] = i;
    }
    for (i=0; i<Ni-1; i++) {
      k = Coefs[i];  m = i;  c = 0;
      for (j=i+1; j<Ni; j++)
  if (Coefs[j] > k) {
    k = Coefs[j];  m = j;  c = 1;
  }
      if (c) {
  Coefs[m] = Coefs[i];  Coefs[i] = k;
  k = InOrder[m];  InOrder[m] = InOrder[i];  InOrder[i] = k;
      }
    }
    for (i++; i<N-1; i++) {
      k = Coefs[i];  m = i;  c = 0;
      for (j=i+1; j<N; j++)
  if (Coefs[j] < k) {
    k = Coefs[j];  m = j;  c = 1;
  }
      if (c) {
  Coefs[m] = Coefs[i];  Coefs[i] = k;
  k = InOrder[m];  InOrder[m] = InOrder[i];  InOrder[i] = k;
      }
    }
  }
}

int gcdmult(m11,m21,pm11)
int m11, m21, *pm11;
{ int m12=1, m22=0, k;

  while (m11 && m21)
    if (m11<m21) {
      k = m21/m11;  m21 -= k*m11;  m22 -= k*m12;
    } else {
      k=m11/m21;  m11 -= k*m21;  m12 -= k*m22;
    }
  if (m11) {
    *pm11 = m12;  return(m11);
  }
  *pm11 = m22;
  return(m21);
}

int alive()
{ int r;  BITS *pd = DeadConfs;

  for (r=0; r<Dead; r++, pd++)
    if ((Mask & *pd) == *pd)  return(0);
  return(1);
}

int aliveineq()
{ int r;  BITS *pd = EDeadCs;

  for (r=0; r<EDead; r++, pd++)
    if ((Mask & *pd) == *pd)  return(0);
  pd = DeadConfs;
  for (r=0; r<Dead; r++, pd++)
    if ((Mask & *pd) == *pd)  return(0);
  return(1);
}

int newval(i0,im,ps,psv,sup)
int i0, im, *ps, *psv, sup;
{ 
  im--;
  if (*psv < sup && ConfVals[im]+1 < CMaxs[im]) {
    ConfVals[im]++;  *ps += Coefs[ConfVars[im]];  (*psv)++;
    return(1);
  }
  if (im) {
    do {
      *psv -= ConfVals[im]-1;  *ps -= (ConfVals[im]-1)*Coefs[ConfVars[im]];
      ConfVals[im--] = 1;
    } while (im >= i0 && (*psv >= sup || ConfVals[im]+1 >= CMaxs[im]));
    if (im < i0)  return(0);
    ConfVals[im]++;  *ps += Coefs[ConfVars[im]];  (*psv)++;
    return(1);
  }
  return(0);
}

int newconfig(n,k,off,ps,psv)
int n, k, off, *ps, *psv;
{ int v = k-1+off;
 
  if (k) {
    if (ConfVars[v] < n-1) {
      Mask &= ~Bit[ConfVars[v]];
      *psv -= ConfVals[v]-1;  *ps -= ConfVals[v]*Coefs[ConfVars[v]];
      ConfVals[v] = 1;  Mask |= Bit[++ConfVars[v]];
      *ps += Coefs[ConfVars[v]];
      return(1);
    }
    if (v > off) {
      do {
  Mask &= ~Bit[ConfVars[v]];  *ps -= Coefs[ConfVars[v--]];
      } while (v>=off && ConfVars[v]==n-k+v-off);
      if (v >= off) {
  Mask &= ~Bit[ConfVars[v]];  *ps -= Coefs[ConfVars[v]];
  Mask |= Bit[++ConfVars[v]];  *ps += Coefs[ConfVars[v]];
  for (v++; v<k+off; v++) {
    Mask |= Bit[(ConfVars[v] = ConfVars[v-1]+1)];
    *ps += Coefs[ConfVars[v]];
  }
  return(1);
      }
    }
  }
  return(0);
}

void slopes()
{  int dz = XDzij, ymax = XMaxi, ddy, dy = XDyij;
   int dt = 1, t = 0, y = ymax-1, k, v;
   
/* this solves equation before equation 4 of \cite{FilgT93 - EPIA}
   NOT the equation after equation 4
*/

   YSlop[ymax] = ymax;  VSlop[ymax] = v = -Eb*ymax;  ZSlop[ymax] = 0;
   SmDy = ymax;   /* the smallest dy */

   if (dy){
     ddy = ymax%dy;     v = -Eb*dy + Ec*dz;   SmDy = dy;
     while (y >= dy){
       YSlop[y] = dy;  VSlop[y] = v;  ZSlop[y--] = dz;  
     }
     while (ddy){
       while (dy-ddy > 0){
   dy -= ddy;  t += dt;
   /* paper notation:  dz = k*dz1, ymax*t <= dy1*k */
   dz = ((t*ymax)/XDyij + 1)*XDzij;   /* overflow! */
   v = -Eb*dy + Ec*dz;   SmDy = dy;
   while (y >= dy){
     YSlop[y] = dy;  VSlop[y] = v;  ZSlop[y--] = dz;
   }
       }
       k = ddy/dy;  ddy = ddy%dy;  dt += k*t;
     }
   }
   VLgDz = v;       
   /* VLgDz is the value corresponding to the m-spacing with the */
   /* largest dz < XMaxj */ 

   v = Ec*XMaxj;
   while (y){
     YSlop[y] = 0;  VSlop[y] = v;  ZSlop[y--] = XMaxj;
   }
 }


void apt(i,j,dz,dy,z2,y1)
int i, j, dz, dy, z2, y1;
{ int z=0, x, k;
  
  while (z+dz < z2) {
    while (y1-dy > 0) {
      y1 -= dy;  z += dz;
      x = Sols[NSols][Eix] = (Eb*y1+Ec*z)/Ea;
      if ((Sols[NSols][i] = y1) == 1)
  if ((Sols[NSols][j] = z) == 1)
    if (x == 1)  DeadConfs[Dead++] = Mask;
    else {
      EDeadCs[EDead++] = Mask;
      LMaxs[Eix] = newel(x,Mask,LMaxs[Eix]);
    }
  else {
    ELMaxs[j] = newel(z,Mask,ELMaxs[j]);
    if (x == 1)  LMaxs[j] = newel(z,Mask,LMaxs[j]);
  }
      else if ((Sols[NSols][j] = z) == 1) {
  ELMaxs[i] = newel(y1,Mask,ELMaxs[i]);
  if (x == 1)  LMaxs[i] = newel(y1,Mask,LMaxs[i]);
      }
      NZSol[NSols++] = Mask;
    } 
    if (!(y1-dy))  return;
    k = dy/y1;  dy = dy%y1;  dz += k*z;
  }
}


int deadfilter(i,j,y,z)
int i, j, y, z;
{ int  ge, c, no = 0, none, *pvr, *pvl, x, *pcs, r;
  BITS *nz, bi = Bit[i], bj = Bit[j];

  for (r=CCSol; r<NSols; r++)
    if (Sols[r][i] <= y && Sols[r][j] <= z) {
      ge = 1;
      for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++)
        if (Sols[r][*pvr++] > *pvl++) {
          ge = 0;
          break;

        }
#ifdef COUNTCAND
      if (ge) {
        RejF++;  return(0);
      }
#else
      if (ge) return(0);
#endif
    }
  for (r=0, pcs=CompSols; r<LastCmpSol; r++, pcs++)
    if (Sols[*pcs][i] <= y && Sols[*pcs][j] <= z) {
      ge = 1;  nz = NZSol+*pcs;
      for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++, pvr++, pvl++)
  if ((*nz&Bit[*pvr]) && Sols[*pcs][*pvr] > *pvl) {
    ge = 0;
    break;
  }
#ifdef COUNTCAND
      if (ge) {
  RejF++;  return(0);
      }
#else
      if (ge) return(0);
#endif
    }
  if (NextCmpSol > NONE) {
    if (NextCmpSol == ALL)  NextCmpSol = EL2Sol;
    for (nz=NZSol+NextCmpSol; NextCmpSol<CLSol; NextCmpSol++, nz++)
      if ((Mask | *nz) == Mask) {
  if (!(*nz&bi))  Sols[NextCmpSol][i] = 0;
  if (!(*nz&bj))  Sols[NextCmpSol][j] = 0;
  *pcs++ = NextCmpSol;  LastCmpSol++;
  if (Sols[NextCmpSol][i] <= y && Sols[NextCmpSol][j] <= z) {
    ge = 1;
    for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++, pvr++, pvl++)
      if ((*nz&Bit[*pvr]) && Sols[NextCmpSol][*pvr] > *pvl) {
        ge = 0;
        break;
      }
    if (ge) {
#ifdef COUNTCAND
      RejF++;  
#endif
      if (NextCmpSol == CLSol)  NextCmpSol = NONE;
      else  NextCmpSol++;
      return(0);
    }
  }
      }
    NextCmpSol = NONE;
  }
  NZSol[NSols] = Mask;
  if ((Sols[NSols][i] = y)==1)  no++;
  if ((Sols[NSols][j] = z)==1)  no++;
  x = z*Ec+y*Eb;
  for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++) {
    x += *pvl*Coefs[*pvr];
    if ((Sols[NSols][*pvr++] = *pvl++)==1)  no++;
    else none = c;
  }
  x /= Ea;
  Sols[NSols++][Eix] = x;
  if (no == ++c)
    if (y > 1)
      if (x == 1)  LMaxs[i] = newel(y,Mask,LMaxs[i]);
      else  ELMaxs[i] = newel(y,Mask,ELMaxs[i]);
    else if (z > 1)
      if (x == 1)  LMaxs[j] = newel(z,Mask,LMaxs[j]);
      else  ELMaxs[j] = newel((CMaxj=z),Mask,ELMaxs[j]);
    else {
      z = CMaxs[none] = ConfVals[none];  j = ConfVars[none];
      if (x == 1)  LMaxs[j] = newel(z,Mask,LMaxs[j]);
      else  ELMaxs[j] = newel(z,Mask,ELMaxs[j]);
    }
  else if (no == ++c) {
    if (x == 1)  DeadConfs[Dead++] = Mask;
    else {
      EDeadCs[EDead++] = Mask;
      LMaxs[Eix] = newel(x,Mask,LMaxs[Eix]);
    }
    return(1);
  }
  return(0);
}

int aptfilter(i,j,v)
int i, j, v;
{ int z, y, dz, dy, maxz = CMaxj, maxy = CMaxi;

  if (maxz > Ea-ESV)  maxz = Ea-ESV;
  if (maxy > Ea-ESV)  maxy = Ea-ESV;
  z = -v/XGCDa2ij*XMultA2ij;  MODASS(z,XDzij);
  y = (-v-z*Ec)*XMultAi/XGCDai;  MODASS(y,XMaxi);
  if (!z) {
    dz = ZSlop[y];
    if (!dz || (z += dz) >= maxz) return(1);
    dy = YSlop[y];  y -= dy;
  } else {
    dy = YSlop[y]; dz = ZSlop[y];
  }
  if (!y) return(1);
  if (!dy){
    if(y < maxy && ESV+z+y <= Ea && deadfilter(i,j,y,z))  return(0);
    return(1);
  }
  do{
    do
      if(y < maxy && ESV+z+y <= Ea && deadfilter(i,j,y,z))  return(0);
    while( (y -= dy) > 0 && (z += dz) < maxz);
    if ( !y || z >= maxz ) return(1);
    y += dy;  dy = YSlop[y];  dz = ZSlop[y];
    y -= dy;  z += dz;
  } while(dy && y && z < maxz);
  return(1);
}


int maxval(i)
int i;
{ LIST l = LMaxs[i];

  i = Max[i];
  while (l != NULL) {
    if ((l->cmask & Mask) == l->cmask && i > l->cmax)
      i = l->cmax;
    l = l->nxt;
  }
  return(i);
}

int maxvaleq(i)
int i;
{ int m = maxval(i);  LIST l = ELMaxs[i];

  while (l != NULL) {
    if ((l->cmask & Mask) == l->cmask && m > l->cmax)  m = l->cmax;
    l = l->nxt;
  }
  if (m > EMax[i])  return(EMax[i]);
  return(m);
}

void ax0bjyj()
{ int i, j, k, n, s, t, v, gb, mb, *pmxy, *pmxz;
  BITS msk, nsk, *bi, *bj;

  Ea = Coefs[Eix];  Mask = msk = Bit[Eix];
  for (i=E0, pmxy=EMax+E0; i<ENj; i++)  *pmxy++ = Ea/GCD[Eix][i];
  EL2Sol = NSols;  EDead = 0;
  for (i=E0, bi=Bit+E0, pmxy=EMax+E0; i<ENj-1; i++, bi++, pmxy++) {
    Eb = Coefs[i];  gb = GCD[Eix][i];  mb = Mult[Eix][i];
    for (j=i+1, bj=Bit+i+1, pmxz=EMax+i+1; j<ENj; j++, pmxz++) {
      Mask = msk|*bi|*bj++;
      if (aliveineq()) {
  Ec = Coefs[j];  t = GCDa2[i][j] = gcdmult(Ec,gb,&(MultA2[i][j]));
  s = Dz[i][j] = gb/t;  v = Ec*mb/t;  MODASS(v,(*pmxy));  Dy[i][j] = v;
  XGCDAj = GCD[Eix][j];  XGCDa2ij = t;
  apt(i,j,s,v,*pmxz,*pmxy);
      }
    }
  }

  for (CVars=1, CLSol=NSols; CVars<=ERhsL-2-(EDead ? 1 : 0);
       CVars++, CLSol=NSols) {
    n = ENj-CVars;
    for (i=E0, bi=Bit+E0; i<n-1; i++) {
      nsk = msk|*bi++;  Eb = Coefs[i];
      XMaxi = EMax[i];  XMultAi = Mult[Eix][i];  XGCDai = GCD[Eix][i];
      for (j=i+1, bj=Bit+i+1; j<n; j++) {
  Mask = nsk|*bj++;
  if (aliveineq()) {
    XMaxj = EMax[j];  Ec = Coefs[j];
    XDzij = Dz[i][j];  XDyij = Dy[i][j];
    XGCDa2ij = GCDa2[i][j];  XMultA2ij = MultA2[i][j];
    XGCDAj = GCD[Eix][j];  XMultAj = Mult[Eix][j];
    s = 0;  ESV = CVars;
    slopes();
    for (k=0, t=j+1; k<CVars; k++, t++) {
      ConfVars[k] = t;  ConfVals[k] = 1;
      s += Coefs[t];  Mask |= Bit[t];
    }
    do
      if (aliveineq()) {
        NextCmpSol = ALL;  LastCmpSol = 0;
        CMaxj = maxvaleq(j);  CMaxi = maxvaleq(i);  CCSol = NSols;
        for (k=0; k<CVars; k++)  CMaxs[k] = maxvaleq(ConfVars[k]);
        /* k must be non-zero here ! */
        do
    if ((v = s%Ea) && !(v % XGCDa2ij))
      k = aptfilter(i,j,v);
        while (k && newval(0,CVars,&s,&ESV,Ea));
      }
    while (newconfig(ENj,CVars,0,&s,&ESV));
  }
      }
    }
  }
}

int deadsfilter(i,j,y,z,x)
int i, j, y, z, x;
{ int  c, *pvr, *pvl, no = 0, none, *pcs, r;
  BITS *nz, bi=Bit[i], bj=Bit[j], bx=Bit[Eix];

  x /= Ea;
  for (r=CCSol; r<NSols; r++)
    if (Sols[r][Eix] <= x && Sols[r][i] <= y && Sols[r][j] <= z) {
      for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++)
        if (Sols[r][*pvr++] > *pvl++)  goto nxtsol;
#ifdef COUNTCAND
      RejF++;  
#endif
      return(0);
    nxtsol:;
    }
  for (r=0, pcs=CompSols; r<LastCmpSol; r++, pcs++)
    if (Sols[*pcs][i] <= y && Sols[*pcs][j] <= z && Sols[*pcs][Eix] <= x) {
      nz = NZSol+*pcs;
      for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++, pvr++, pvl++)
  if ((*nz&Bit[*pvr]) && Sols[*pcs][*pvr] > *pvl)  goto nxtcsol;
#ifdef COUNTCAND
      RejF++; 
#endif
      return(0);
    nxtcsol:;
    }
  if (NextCmpSol > NONE) {
    if (NextCmpSol == ALL)  NextCmpSol = 0;
    for (nz=NZSol+NextCmpSol; NextCmpSol<CLSol; NextCmpSol++, nz++)
      if ((Mask | *nz) == Mask) {
  if (!(*nz&bx))  Sols[NextCmpSol][Eix] = 0;
  if (!(*nz&bi))  Sols[NextCmpSol][i] = 0;
  if (!(*nz&bj))  Sols[NextCmpSol][j] = 0;
  *pcs++ = NextCmpSol;  LastCmpSol++;
  if (Sols[NextCmpSol][i] <= y && Sols[NextCmpSol][j] <= z &&
      Sols[NextCmpSol][Eix] <= x) {
    for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++, pvr++, pvl++)
      if ((*nz&Bit[*pvr]) && Sols[NextCmpSol][*pvr] > *pvl)
        goto nxtssol;
#ifdef COUNTCAND
    RejF++;  
#endif
    if (NextCmpSol == CLSol)  NextCmpSol = NONE;
    else  NextCmpSol++;
    return(0);
  }
      nxtssol:;
      }
    NextCmpSol = NONE;
  }
  if ((Sols[NSols][i] = y) == 1)  no++;
  if ((Sols[NSols][j] = z) == 1)  no++;
  for (c=0, pvr=ConfVars, pvl=ConfVals; c<CVars; c++)
    if ((Sols[NSols][*pvr++] = *pvl++) == 1)  no++;
    else  none = c;
  if ((Sols[NSols][Eix] = x) == 1)  no++;
  NZSol[NSols++] = Mask;
  if (no == CVars+3) {
    DeadConfs[Dead++] = Mask;
    return(1);
  } else if (no == CVars+2)
    if (x > 1)  LMaxs[Eix] = newel(x,Mask,LMaxs[Eix]);
    else if (y > 1)  LMaxs[i] = newel(y,Mask,LMaxs[i]);
    else if (z > 1)  LMaxs[j] = newel(z,Mask,LMaxs[j]);
    else {
      z = CMaxs[none] = ConfVals[none];
      j = ConfVars[none];
      LMaxs[j] = newel(z,Mask,LMaxs[j]);
    }
  return(0);
}


int saptfilter(i,j,vS0)
int i, j, vS0;
{ int y0, z0, czmax = Ec*XMaxj, dymin = SmDy, k, maxy = CMaxi, 
    maxz = CMaxj, bymax = Eb*XMaxi, zplus = (-vS0-1)/Ec+1;

  /* upper bounds for y and z */

  if (maxz > LhsSup-RhsSV)  maxz = LhsSup-RhsSV;
  if (maxy > LhsSup-RhsSV)  maxy = LhsSup-RhsSV;

  if ( vS0 <= (1-maxz)*Ec+(1-maxy)*Eb) return(1);  /* every sol is negative */

  /* First solution */

  z0 = -vS0/XGCDa2ij*XMultA2ij;  MODASS(z0,XDzij);
  y0 = (-vS0-z0*Ec)*XMultAi/XGCDai;  MODASS(y0,XMaxi);

  if (!z0) 
    if (vS0 + Eb*y0 + Ec*z0 >= 0) { 
      if (maxy > y0) maxy = y0;
    } else {
#ifdef COUNTCAND
      Rej0++;            /* not counted in posaptfilter */
#endif
      z0 += XDzij; 
      if ((y0 -= XDyij) < 0) y0 += XMaxi;
    }
  
  
  if ((vS0 += (Eb*y0 + Ec*z0)) < 0) {
    k =  (-vS0-1)/bymax + 1;    /* spacing for positiveness */ /* ceil */
    vS0 += (bymax*k);  y0 += (XMaxi*k);
#ifdef COUNTCAND
    Rej0++;
#endif
  } 
  /* It's possible that y0 >= maxy. The algorithm to ensure y0 < maxy */
  /* would much resemble the previous saptfilter version */ 

  if (zplus > maxz) zplus = maxz;

  /* Other solutions */ 

  while (z0 < zplus && y0 > dymin) {
    if (y0 < maxy && vS0 && RhsSV+z0+y0<=LhsSup)
      deadsfilter(i,j,y0,z0,vS0);
#ifdef COUNTCAND
    else Rej0++;
#endif
    if (VLgDz + vS0 < 0) {   
      /* no m-spacing but (0,XMaxj) applies; all of them are negative */
      /* and greater in absolute value than vS0 */ 
      k = (-VLgDz-vS0-1)/czmax + 1;
      if ( (z0 += k*XMaxj) >= maxz ) {
#ifdef COUNTCAND
  Rej0++;
#endif
  return (1);
      }
      vS0 += k*czmax;   
    }

    for (k = (y0>=XMaxi? XMaxi : y0); VSlop[k] + vS0 < 0; k = YSlop[k-1]) ;   
    vS0 += VSlop[k];  y0 -= YSlop[k];  z0 += ZSlop[k];  
    if (y0 >= XMaxi && vS0 >= bymax) {
      if ((k = vS0/bymax) > y0/XMaxi) k = y0/XMaxi;
      vS0 -= k*bymax;   y0 -= k*XMaxi;
    }
  }

  while (z0 < maxz && y0 > dymin) {     /* non-negative values */
    if (y0 < maxy && vS0 && RhsSV+z0+y0<=LhsSup)
      deadsfilter(i,j,y0,z0,vS0);
#ifdef COUNTCAND
    else Rej0++;
#endif
    vS0 += VSlop[y0];  z0 += ZSlop[y0];  y0 -= YSlop[y0];
  }

  if (z0 < maxz && y0 && y0 < maxy && vS0 && RhsSV+z0+y0<=LhsSup) {
    if (deadsfilter(i,j,y0,z0,vS0)) return (0);
  }
#ifdef COUNTCAND
   else  Rej0++;
#endif
  return (1);
}

int posaptfilter(i,j,r)
int i, j, r;
{ int z, y, maxz = CMaxj, maxy = CMaxi, dz, dy;
  /* not the same as in slopesV5i.c because the former m-spacing (0,0)
     has been replaced by (0,XMaxj) in slopes */

  if (maxz > LhsSup-RhsSV)  maxz = LhsSup-RhsSV;
  if (maxy > LhsSup-RhsSV)  maxy = LhsSup-RhsSV;
  z = -r/XGCDa2ij*XMultA2ij;  MODASS(z,XDzij);
  y = (-r-z*Ec)*XMultAi/XGCDai;  MODASS(y,XMaxi);
  if (!z){
    if ( !(dy = YSlop[y]) || (z += (dz = ZSlop[y])) > maxz) return(1);
    y -= dy;
  }else { 
    dy = YSlop[y]; dz = ZSlop[y];
  }
  if (!y) return(1);
  if (!dy) {
    if (y < maxy && RhsSV+z+y <= LhsSup) {
      if (deadsfilter(i,j,y,z,Eb*y+Ec*z+r)) return(0);
    }
#ifdef COUNTCAND
      else  Rej0++;
#endif
    return(1);
  }
  do{
    do {
      if (y < maxy && RhsSV+z+y <= LhsSup) {
  if (deadsfilter(i,j,y,z,Eb*y+Ec*z+r)) return(0);
      }
#ifdef COUNTCAND
      else  Rej0++;
#endif
    } while ((y -= dy) > 0 && (z += dz) < maxz);
    if ( z >= maxz || !y) return(1);
    y += dy;  dy = YSlop[y];  dz = ZSlop[y];
    y -= dy;  z += dz;
  } while (z < maxz && dy && y);
  return(1);
}

void axbjyj()
{ int i, j, k, n, s, sum, t, *pmxy, *pmxz, rhscv, maxbc;
  BITS msk = Mask, nsk, *bi, *bj;

  CVars = LhsL;  sum = LhsSum;
  for (i=E0, pmxy=EMax+E0; i<ENj; i++)  *pmxy++ = Ea/GCD[Eix][i];
  for (rhscv=0, CLSol=NSols; rhscv<=ERhsL-2 && CVars+3<=N+(Dead ? 1 : 0);
       rhscv++, CVars++, CLSol=NSols) {
    n = ENj-rhscv;
    for (i=E0, bi=Bit+E0, pmxy=EMax+E0; i<n-1; i++, bi++, pmxy++) {
      nsk = msk|*bi;
      Eb = Coefs[i];  XGCDai = GCD[Eix][i];  XMultAi = Mult[Eix][i];
      XMaxi = *pmxy;
      for (j=i+1, bj=Bit+i+1, pmxz=EMax+i+1; j<n; j++, pmxz++) {
  Mask = nsk|*bj++;
  if (alive()) {
    if ((Ec = Coefs[j]) > Eb)  maxbc = Ec;
    else  maxbc = Eb;
    XMaxj = *pmxz;  XGCDa2ij = gcdmult(Ec,XGCDai,&XMultA2ij);
    XDzij = XGCDai/XGCDa2ij;
    XDyij = Ec*XMultAi/XGCDa2ij;  MODASS(XDyij,XMaxi);
    XGCDAj = GCD[Eix][j];  XMultAj = Mult[Eix][j];
    s = 0;  RhsSV = rhscv;
    slopes();
    for (k=LhsL, t=j+1; k<CVars; k++, t++) {
      ConfVars[k] = t;  ConfVals[k] = 1;
      s += Coefs[t];  Mask |= Bit[t];
    }
    do
      if (alive()) {
        NextCmpSol = ALL;  LastCmpSol = 0;
        CCSol = NSols;  RhsSup = maxbc;
        CMaxi = maxval(i);  CMaxj = maxval(j);
        for (k=0; k<LhsL; k++)  CMaxs[k] = maxval(ConfVars[k]);
        for (; k<CVars; k++) {
    CMaxs[k] = maxval(ConfVars[k]);
    if (Coefs[ConfVars[k]] > RhsSup)  RhsSup = Coefs[ConfVars[k]];
        }
        /* k not zero here */
        do
    if (! ((t = s-LhsSum)%XGCDa2ij))
      if (t < 0)  k = saptfilter(i,j,t);
    else if (t)  k = posaptfilter(i,j,t);
        while (k &&
           ((rhscv && newval(LhsL,CVars,&s,&RhsSV,LhsSup)) ||
      newval(0,LhsL,&LhsSum,&LhsSV,RhsSup)));
        for (t=0; t<LhsL; t++)  ConfVals[t] = 1;
        LhsSum = sum;  LhsSV = LhsL;
      }
    while (newconfig(ENj,rhscv,LhsL,&s,&RhsSV));
  }
      }
    }
  }
}

int main(argc,argv)
int argc;
char *argv[];
{ int t, mi, mj, mm, *pi, *pj, *pmi, i, j, k, many = 0, sols = 0, times = 1;
  BITS *nz, msk, *pb, *pbj;
  LIST *pl, *pel;
  double time; 

  while (--argc)
    switch (*argv[argc]) {
    case 'm':  many = 1;  break;
    case 's':  sols = 1;  break;
#ifndef COUNT
    case 't':  times = 100;  break;
    default:  printf("s: print solutions, m: solve many problems, t: *100\n");
#else
    default:  printf("s: print solutions, m: solve many problems\n");
#endif
      exit(1);
    }
  do {
    inputdata();
    N = Ni+Nj;  OutOrder = 0;
    if (Nj < Ni) {
      for (i=0, pi=Coefs, pj=Coefs+Ni; i<Ni; i++) {
  t = *pi;  *pi++ = *pj;  *pj++ = t;
      }
      for (i=0, pi=Coefs+Nj, pj=Coefs+Ni; i<Ni; i++)  *pi++ = *pj++;
      t = Ni;  Ni = Nj;  Nj = t;
    }
    for (t=0; t<times; t++) {
      LhsSup = RhsSup = 0;
      for (i=0, pi=Coefs; i<Ni; i++, pi++)
  if (*pi > LhsSup)  LhsSup = *pi;
      for (; i<N; i++, pi++)
  if (*pi > RhsSup)  RhsSup = *pi;
      ordercoeffs();
      NSols = Dead = 0;
#ifdef COUNTCAND
      Rej0 = RejF = 0;
#endif
      pb = Bit+Ni;  pl = LMaxs+Ni;  pel = ELMaxs+Ni;  pmi = Max+Ni;
      msk = (1<<Ni);
      for (j=0; j<Nj; j++, msk=(msk<<1)) {
  *pb++ = msk;  *pl++ = *pel++ = NULL;  *pmi++ = LhsSup;
      }
      pb = Bit;  pl = LMaxs;  pel = ELMaxs;  pmi = Max;
      msk = 1;  nz = NZSol;  pi = Coefs;
      for (i=0; i<Ni; i++, pi++, pb++, msk=(msk<<1)) {
  *pb = msk;  *pl++ = *pel++ = NULL;  *pmi++ = RhsSup;
  for (j=Ni, pj=Coefs+Ni, pbj=Bit+Ni; j<N; j++, nz++) {
    *nz = *pb | *pbj++;
    k = GCD[i][j] = GCD[j][i] = gcdmult(*pj,*pi,&mm);
    Mult[i][j] = mm;  Mult[j][i] = (k-mm**pj)/(*pi);
    mi = *pi/k;  mm = *pj++/k;
    if ((Sols[NSols][i] = mm) == 1) {
      if ((Sols[NSols][j] = mi) == 1) {
        DeadConfs[Dead++] = *nz;
      } else {
        LMaxs[j] = newel(mi,*nz,LMaxs[j]);
      }
    } else if ((Sols[NSols][j] = mi) == 1) {
      LMaxs[i] = newel(mm,*nz,LMaxs[i]);
    }
    NSols++;
  }
      }
      L2Sol = NSols;
      if (Ni > 1) {
  E0 = 0;  ENj = Ni;  ERhsL = Ni;
  for (Eix=Ni; Eix<N; Eix++)  ax0bjyj();
      }
      E0 = Ni;  ENj = N;  ERhsL = Nj;
      for (Eix=0; Eix<Ni; Eix++)  ax0bjyj();
      for (LhsL=1; LhsL<Ni; LhsL++)  /* LhsL: lhs cvars */
  for (Eix=0, pb=Bit; Eix<Ni-LhsL; Eix++) {
    LhsSum = 0;  Mask = *pb++;  LhsSV = LhsL;  Ea = Coefs[Eix];
    for (mj=0, j=Eix+1, pbj=Bit+Eix+1; mj<LhsL; mj++, j++) {
      ConfVars[mj] = j;  ConfVals[mj] = 1;
      LhsSum += Coefs[j];  Mask |= *pbj++;
    }
    do {
      LhsSup = Ea;
      for (mj=0; mj<LhsL; mj++)
        if (Coefs[ConfVars[mj]] > LhsSup)  LhsSup = Coefs[ConfVars[mj]];
      msk = Mask;  axbjyj();  Mask = msk;
    } while (newconfig(Ni,LhsL,0,&LhsSum,&LhsSV));
  }
    }
    /* stop clock */
    if (sols)
      for (k=0; k<NSols; k++)  printsol(k);
    if (N) {
      //      printf("Problem:\t");
      //     for (k=0; k<Ni; k++)  printf(" %d",Coefs[actual(k)]);
      //      printf(" (=) ");
      //      for (k=Ni; k<N; k++)  printf(" %d",Coefs[actual(k)]);
      printf("%d\n",NSols);
#ifdef COUNTCAND
      //      printf("\n%s:\tNo. sols = %d\n",Version,NSols);
      //      printf("\tNo. lim. rejected = %f, No. filtered = %f\n",Rej0,RejF);
#else
      //      printf("\n%s:\tNo. sols = %d,  cputime = %f/%d sec = %f\n",
      //       Version,NSols,time,times,time/times);
#endif
      fflush(stdout);
    }
  } while (many && N);
  return(0);
}
    
