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

#include <array>
#include <future>
#include <thread>
#include <cassert>
#include <NTL/lzz_pXFactoring.h>
#include <SDL2/SDL.h>
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "matrix.h"

template<class type> 
class TransformMatrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< RX > > data;

public:
  ~TransformMatrix() { }

  TransformMatrix(const EncryptedArray& _ea) : ea(_ea) {
    long n = ea.size();
    long d = ea.getDegree();

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.resize(n);
    for (long i = 0; i < n; i++) {
      data[i].resize(n);
      for (long j = 0; j < n; j++) {
        bool color_row = j == n - 1;

        if (color_row) {
          clear(data[i][j]);
        } else {
          data[i][j] = 1;
        }
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
buildTransformMatrix(const EncryptedArray& ea)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new TransformMatrix<PA_GF2>(ea);
    }

    case PA_zz_p_tag: {
      return new TransformMatrix<PA_zz_p>(ea);
    }

    default: return 0;
  }
}

typedef std::array<float, 4> Vec4;

struct EncodeState
{
  Vec<ZZX>& v;
  shared_ptr<EncryptedArray> ea;
  std::vector<GF2X> slots;
  int slot;
  int index;
  
  EncodeState(shared_ptr<EncryptedArray> ea, Vec<ZZX>& v)
  : ea(ea), v(v), slot(0), index(0) {
        assert(ea->size() >= 4);
        assert(ea->size() % 4 == 0);
        slots = std::vector<GF2X>(ea->size(), GF2X::zero());
  }

  void EncodeVector(const Vec4& vec)
  {
    const EncryptedArrayDerived<PA_GF2>& ea2 = ea->getDerived(PA_GF2());
    
    if (slot + 4 >= ea2.size()) {
        slot = 0;
        ea2.encode(v[index], slots);
        slots = std::vector<GF2X>(ea2.size(), GF2X::zero());
        ++index;
    }

#if 1
    cout << "{";
#endif

    // Convert float (-1.0 to 1.0) to 16bit signed magnitude int.
    for (int i = 0; i < 4; ++i)
    {
        assert(vec[i] <= 1.0 && vec[i] >= -1.0);
        float fval = vec[i];
        int16_t ival = 0;
        double fval_range = abs(fval) * (double)0x7FFF;
        ival = round(fval_range);
        if (fval < 0.0) {
            ival += 0x8000;
        }
        GF2XFromBytes(slots[slot++], (const unsigned char*)&ival, 2);
#if 1
        int16_t dbg_ival = 0;
        cout << fval << "->";
        cout << ival << "=";
        BytesFromGF2X((unsigned char*)&dbg_ival, slots[slot - 1], 2);
        cout << dbg_ival << ",";
        assert(dbg_ival == ival);
#endif
    }
#if 1
    cout << "}\n";
#endif
  }
};

void GetColor(const Vec4& vertex, Vec4& out)
{
  out[0] = vertex[0] * 255;
  out[1] = vertex[1] * 255;
  out[2] = vertex[2] * 255;
  out[3] = vertex[3] * 255;
}

void DecodeVectors(shared_ptr<EncryptedArray>& ea, ZZX& encoded, vector<Vec4>& out)
{
  const EncryptedArrayDerived<PA_GF2>& ea2 = ea->getDerived(PA_GF2());
  
  std::vector<GF2X> slots;
  ea2.decode(slots, encoded);
 
  for (int j = 0; j < ea2.size() / 4; ++j) {
#if 1
        cout << "{";
#endif

    // Convert 16bit signed magnitude int to float (-1.0 to 1.0).
    for (int i = 0; i < 4; ++i)
    {
        Vec4 tmp;
        int16_t ival = 0;
        BytesFromGF2X((unsigned char*)&ival, slots[(j * 4) + i], 2);
        
        if (ival & 31) {
            ival -= 0x8000;
            tmp[i] = -1.0;
        } else {
            tmp[i] = 1.0;
        }
      
        tmp[i] *= ival / 0x7FFF;
#if 1
        cout << tmp[i] << "<-";
        cout << ival << ",";
#endif
        out.push_back(tmp);
    }
#if 1
        cout << "}\n";
#endif
  }
}

struct State
{
FHEcontext context;
shared_ptr<FHESecKey> secretKey;
shared_ptr<FHEPubKey> publicKey;
shared_ptr<EncryptedArray> ea;

Vec<ZZX> vertexdata;
std::vector<Ctxt> ctxt_vertexdata;

std::vector<Ctxt> ctxt_framebuffer;
std::vector<Ctxt> ctxt_clearcolor;

int vertex_stride;
int triangle_stride;
int num_triangles;
int width;
int width_pixelset;
int pixelset_size;
int height;

shared_ptr<PlaintextMatrixBaseInterface> transform;

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

  secretKey = shared_ptr<FHESecKey>(new FHESecKey(context));
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
  ea = shared_ptr<EncryptedArray>(new EncryptedArray(context, G));
  cout << "done\n";
  cout << ea->size() << "\n";
  cout << ea->getDegree()  << "\n";
 
  // A triangle
  cout << "Input triangle\n";
  vertex_stride = 4;
  triangle_stride = 4;
  num_triangles = 1;
  
  vertexdata.SetLength(1); // Calculate this
  
  {
    EncodeState tmp(ea, vertexdata);
    tmp.EncodeVector({-0.37234,-1.0,0.0,0.0}); // Bottom Left
    tmp.EncodeVector({0.0,1.0,0.0,0.0}); // Top Middle
    tmp.EncodeVector({1.0,-1.0,0.0,0.0}); // Bottom Right
    tmp.EncodeVector({1.0,0.0,1.0,1.0}); // Per triangle color -- change to per vertex
  
    // encrypt the trangle verticies
    for (int i = 0; i < vertexdata.length(); ++i) {
      ctxt_vertexdata.push_back(Ctxt(*publicKey));
      publicKey->Encrypt(ctxt_vertexdata[i], vertexdata[i]);
    }
   
    // test decryption of the trangle verticies
    for (int i = 0; i < vertexdata.length(); ++i) {
      ZZX v1;
      secretKey->Decrypt(v1, ctxt_vertexdata[i]);

      vector<Vec4> a;
      vector<Vec4> b;
      DecodeVectors(ea, v1, a);
      DecodeVectors(ea, vertexdata[i], b);
      assert(a.size() == b.size());
      for (int j = 0; j < a.size(); ++j) {
          if (a[j] != b[j]) {
             cout << "Fail!!\n";
          }
      }
    }
  }
  
  // Constants just clear color for now
  cout << "Input Constants\n";
  Vec<ZZX> constants;
  constants.SetLength(1);
  EncodeState tmp(ea, constants);
  tmp.EncodeVector({0.0,0.0,0.0,1.0});
  ctxt_clearcolor.push_back(Ctxt(*publicKey));
  publicKey->Encrypt(ctxt_clearcolor[0], constants[0]);
  
  height = 64;
  width = 64;
  pixelset_size = ea->size() / 4;
  width_pixelset = width / pixelset_size;
  
  for (long i = 0; i < height * width_pixelset; ++i) {
    ctxt_framebuffer.push_back(ctxt_clearcolor[0]);
  }
    
  // vertex shader rotate matrix
  //transform = shared_ptr<PlaintextMatrixBaseInterface>(buildTransformMatrix(*ea));
}

};

void Render(State& state)
{
    cout << "Render Scene\n";

    // Clear framebuffer
    for (int y = 0; y < state.height; ++y) {
      for (int x_ = 0; x_ < state.width_pixelset; ++x_) {
         const int ps = (y * state.width_pixelset) + x_;
         state.ctxt_framebuffer[ps] = state.ctxt_clearcolor[0];
      }
    }
    
    /*for (int i = 0; i < state.num_triangles; ++i) {
      // Begin Vertex Shader {
      //   mat_mul(*state.ea, state.ctxt_vertexdata[1], *transform);  // transform the triangles
      Ctxt color_out = state.ctxt_vertexdata[(i * 4) + 3];
      // } End Vertex Shader
      
        
      // For each pixel fill/rasterize triangle
      for (int y = 0; y < state.height; ++y) {
        for (int x_ = 0; x_ < state.width_pixelset; ++x_) {
          const int ps = (y * state.width_pixelset) + x_;
          //cout << "PixelSet " << ps << "\n";
          
          // Begin Fragment Shader {
          Ctxt color_in = color_out;
          
          state.ctxt_framebuffer[ps] = color_in;  // blend current pixel with new
          // } End Fragment Shader
        }
      }
    }*/
    
    cout << "\nDone\n\n";
}

void CopyToFramebuffer(State& state, SDL_Renderer *renderer)
{
    cout << "Copy To Framebuffer\n";
    
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    
    for (int y = 0; y < state.height; ++y) {
      for (int x_ = 0; x_ < state.width_pixelset; ++x_) {
        // Decrypt framebuffer pixelsets
        ZZX pixelset;
        const int ps = (y * state.width_pixelset) + x_;
        //cout << "PixelSet " << ps << "\n";
        state.secretKey->Decrypt(pixelset, state.ctxt_framebuffer[ps]);
          
        vector<Vec4> vectors;
        DecodeVectors(state.ea, pixelset, vectors);
        
        // For each pixel in pixel "copy" to screen
        for (int _x = 0; _x < state.pixelset_size; ++_x) {
          const int x = (x_ * state.pixelset_size) + _x;
          
          Vec4 color;
          GetColor(vectors[_x], color);
          
          SDL_SetRenderDrawColor(renderer, color[0],
                                           color[1],
                                           color[2], 255);
            
          SDL_RenderDrawPoint(renderer, x, y);
        }
      }
    }
    
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
  argmap["k"] = "32";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long k = atoi(argmap["k"]);
 
  setTimersOn();
  long m = FindM(k, 4, 3, 2, 16, 4, 771, true);
  State state(m, 2, 1, 16, 4);

  int done;
  SDL_Window *window;
  SDL_Renderer *renderer;

  if (SDL_Init(SDL_INIT_VIDEO) < 0) {
    return -1;
  }

  SDL_CreateWindowAndRenderer(64, 64, 0, &window, &renderer);
  if (window == NULL || renderer == NULL) {
    return -1;
  }

  bool running = true;
  while (running) {
    SDL_Event event;
    if (SDL_PollEvent(&event)) {
      if (event.type == SDL_QUIT) {
         running = false;
      }
    
      if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_SPACE) {
        }
      }
    }
    
    std::future<void> render = std::async(std::launch::async, [&state](){ Render(state); });
    std::chrono::milliseconds timeout(100);
    while (render.wait_for(timeout) == std::future_status::timeout) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
            abort(); // Hard quit for now.
        }
    }
    
    std::future<void> blit = std::async(std::launch::async, [&state, &renderer](){ CopyToFramebuffer(state, renderer); });
    while (blit.wait_for(timeout) == std::future_status::timeout) {
        if (SDL_PollEvent(&event) && event.type == SDL_QUIT) {
            abort(); // Hard quit for now.
        }
    }
    
    SDL_RenderPresent(renderer);
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
