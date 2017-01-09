/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/* Test_matmul.cpp - Testing the functionality of multiplying an encrypted
 * vector by a plaintext matrix, either over the extension- or the
 * base-field/ring.
 */

#include <cassert>
#include <NTL/lzz_pXFactoring.h>
#include <SDL2/SDL.h>
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "matrix.h"

template<class type> 
class RandomMatrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< RX > > data;

public:
  ~RandomMatrix() { /*cout << "destructor: random matrix\n";*/ }

  RandomMatrix(const EncryptedArray& _ea) : ea(_ea) { 
    long n = ea.size();
    long d = ea.getDegree();

    long bnd = 2*n; // non-zero with probability 1/bnd

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.resize(n);
    for (long i = 0; i < n; i++) {
      data[i].resize(n);
      for (long j = 0; j < n; j++) {
        bool zEntry = (RandomBnd(bnd) > 0);

        if (zEntry)
          clear(data[i][j]);
        else
          random(data[i][j], d);
      }
    }
  }

  virtual const EncryptedArray& getEA() const {
    return ea;
  }

  virtual bool get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};

PlaintextMatrixBaseInterface *
buildRandomMatrix(const EncryptedArray& ea)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix<PA_GF2>(ea);
    }

    case PA_zz_p_tag: {
      return new RandomMatrix<PA_zz_p>(ea);
    }

    default: return 0;
  }
}

typedef float Vec4[4];

void AddVertex(std::vector<long>& buffer, const Vec4& vec)
{
  // For now convert float (-1.0 to 1.0) to 32bit signed magnitude int.
  for (int i = 0; i < 4; ++i)
  {
      assert(vec[i] <= 1.0 && vec[i] >= -1.0);
      float fval = vec[i];
      int32_t ival = 0;
      double fval_range = abs(fval) * (double)0x7FFFFFFF;
      ival = round(fval_range);
      if (fval < 0.0) {
          ival += 0x80000000;
      }
      
      for (int i = 0; i < 32; ++i)
      {
          buffer.push_back(ival & 1);
          cout << (ival & 1);
          ival >>= 1;
      }
      cout << "\n";
  }
  cout << "\n";
}

struct State
{
FHEcontext context;
FHESecKey *secretKey;
FHEPubKey *publicKey;
EncryptedArray *ea;

std::vector<NewPlaintextArray> v;
std::vector<Ctxt> ctxt;

std::vector<Ctxt> ctxt_framebuffer;
std::vector<Ctxt> ctxt_clearcolor;

int num_triangles;
int width_pixelset;
int height;

State(long m, long p, long r, long d, long L)
: context(m, p, r)
{
  cout << "*** TestIt: m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", L=" << L
       << endl;

  buildModChain(context, L, /*c=*/3);

  context.zMStar.printout();
  cout << endl;

  secretKey = new FHESecKey(context);
  publicKey = secretKey;
  secretKey->GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  cout << "G = " << G << "\n";
  cout << "generating key-switching matrices... ";
  addSome1DMatrices(*secretKey); // compute key-switching matrices that we need
  addFrbMatrices(*secretKey); // compute key-switching matrices that we need
  cout << "done\n";

  cout << "computing masks and tables for rotation...";
  ea = new EncryptedArray(context, G);
  cout << "done\n";
  
  // A triangle
  cout << "Input triangle\n";
  v.push_back(NewPlaintextArray(*ea));
  std::vector<long> triangle;
  AddVertex(triangle, {-0.37234,-1.0,0.0,0.0}); // Bottom Left
  AddVertex(triangle, {0.0,1.0,0.0,0.0}); // Top Middle
  AddVertex(triangle, {1.0,-1.0,0.0,0.0}); // Bottom Right
  AddVertex(triangle, {1.0,0.0,1.0,1.0});
  triangle.resize(ea->size());
  encode(*ea, v[0], triangle);
  num_triangles = 1;
  
  // encrypt the trangle vector
  ctxt.push_back(Ctxt(*publicKey));
  ea->encrypt(ctxt[0], *publicKey, v[0]);
  
  // test decryption of the trangle vector
  NewPlaintextArray v1(*ea);
  ea->decrypt(ctxt[0], *secretKey, v1);

  if (!equals(*ea, v[0], v1)) {
    cout << "Fail!!\n";
  }
  
  NewPlaintextArray clearcolor(*ea);
  ctxt_clearcolor.push_back(Ctxt(*publicKey));
  ea->encrypt(ctxt_clearcolor[0], *publicKey, clearcolor);
  
  height = 256;
  width_pixelset = 16; // 16 slots per ctxt (pixelset) = 256 pixels row.
  
  for (long i = 0; i < height * width_pixelset; ++i) {
    ctxt_framebuffer.push_back(ctxt_clearcolor[0]);
  }
}

};

void Render(State& state)
{
    // plaintext rotation matrix
    shared_ptr<PlaintextMatrixBaseInterface> ptr(buildRandomMatrix(*state.ea));

    // Clear framebuffer
    for (int y = 0; y < state.height; ++y) {
      for (int x_ = 0; x_ < state.width_pixelset; ++x_) {
        const int i = (y * state.width_pixelset) + x_;
         state.ctxt_framebuffer[i] = state.ctxt_clearcolor[0];
      }
    }
    
    for (int i = 0; i < state.num_triangles; ++i) {
      //mat_mul(*state.ea, state.ctxt[1], *ptr);  // rotate the triangles
      
      // For each pixel fill/rasterize triangle
      for (int y = 0; y < state.height; ++y) {
        for (int x_ = 0; x_ < state.width_pixelset; ++x_) {
        }
      }
    }
}

void CopyToFramebuffer(State& state, SDL_Renderer *renderer)
{
    cout << "Copy To Framebuffer\n";
    
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    
    for (int y = 0; y < state.height; ++y) {
      for (int x_ = 0; x_ < state.width_pixelset; ++x_) {
        // Decrypt framebuffer pixel - currently there are wasted slots
        NewPlaintextArray pa(*state.ea);
        const int i = (y * state.width_pixelset) + x_;
        state.ea->decrypt(state.ctxt_framebuffer[i], *state.secretKey, pa);
        
        std::vector<long> pixelset;
        decode(*state.ea, pixelset, pa);
    
        for (int _x = 0; _x < state.width_pixelset; ++_x) {
          const int x = (x_ * state.width_pixelset) + _x;
          const int color = pixelset[_x];
          
          SDL_SetRenderDrawColor(renderer, (color & 0xFF000000) >> 24,
                                           (color & 0x00FF0000) >> 16,
                                           (color & 0x0000FF00) >> 8, 255);
          SDL_RenderDrawPoint(renderer, x, y);
        }
      }
    }
    SDL_RenderPresent(renderer);
    
    cout << "\nDone\n\n";
}

void usage(char *prog) 
{
  cout << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cout << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cout << "  e.g, 'm=2047 p=2 L=4'\n\n";
  cout << "  m defines the cyclotomic polynomial Phi_m(X)\n";
  cout << "  p is the plaintext base [default=2]" << endl;
  cout << "  r is the lifting [default=1]" << endl;
  cout << "  d is the degree of the field extension [default==1]\n";
  cout << "    (d == 0 => factors[0] defined the extension)\n";
  cout << "  L is the # of primes in the modulus chain [default=4]\n";
  exit(0);
}

/* Testing the functionality of multiplying an encrypted vector by a plaintext
 * matrix, either over the extension- or the base-field/ring.
 */
int main(int argc, char *argv[])
{
  argmap_t argmap;
  argmap["m"] = "8191";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["d"] = "1";
  argmap["L"] = "4";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long m = atoi(argmap["m"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long L = atoi(argmap["L"]);

  setTimersOn();
  State state(m, p, r, d, L);

  int done;
  SDL_Window *window;
  SDL_Renderer *renderer;

  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    return -1;
  }

  SDL_CreateWindowAndRenderer(256, 256, 0, &window, &renderer);
  if (window == NULL || renderer == NULL) {
    return -1;
  }

  bool running = true;
  while (running) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT) {
         running = false;
      }
    
      if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_SPACE) {
        }
      }
    }
    
    Render(state);
  
    CopyToFramebuffer(state, renderer);
  }

  //SDL_DestroyTexture(texture);
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  SDL_Quit();
  
  cout << endl;
  printAllTimers();
  cout << endl;
  
  return 0;
}
