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

typedef std::array<int16_t, 4> Vec4;
#define DEBUG_ENCODE 1
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

  void FinishSlots()
  {
    const EncryptedArrayDerived<PA_GF2>& ea2 = ea->getDerived(PA_GF2());
    
    slot = 0;
    slots = std::vector<GF2X>(ea2.size(), GF2X::zero());
    ++index;
  }
  
  void EncodeVector(const Vec4& vec)
  {
    const EncryptedArrayDerived<PA_GF2>& ea2 = ea->getDerived(PA_GF2());
    
    if (slot + 4 >= ea2.size()) {
        FinishSlots();
    }

#if DEBUG_ENCODE
    cout << "{";
    std::array<int16_t, 4> dbg_ivals;
#endif

    for (int i = 0; i < 4; ++i)
    {
        int16_t ival = vec[i];
        GF2XFromBytes(slots[slot++], (const unsigned char*)&ival, 2);
#if DEBUG_ENCODE
        int16_t dbg_ival = 0;
        cout << ival << "=";
        BytesFromGF2X((unsigned char*)&dbg_ival, slots[slot - 1], 2);
        cout << dbg_ival << ",";
        assert(dbg_ival == ival);
        dbg_ivals[i] = dbg_ival;
#endif
    }
    
    ea2.encode(v[index], slots);
      
#if DEBUG_ENCODE
    cout << "}\n";
    
    std::vector<GF2X> tmp;
    ea2.decode(tmp, v[index]);
    
    for (int i = 0; i < 4; ++i)
    {
        int16_t dbg_ival = 0;
        BytesFromGF2X((unsigned char*)&dbg_ival, tmp[i], 2);
        assert(dbg_ival == dbg_ivals[i]);
        cout << dbg_ival << "=";
        cout << dbg_ivals[i] << ",";
    }
    cout << "\n";
#endif
  }
};

void GetColor(const Vec4& vertex, Vec4& out)
{
  out[0] = round(vertex[0] / 255);
  out[1] = round(vertex[1] / 255);
  out[2] = round(vertex[2] / 255);
  out[3] = round(vertex[3] / 255);
}

void DecodeVectors(shared_ptr<EncryptedArray>& ea, ZZX& encoded, vector<Vec4>& out)
{
  const EncryptedArrayDerived<PA_GF2>& ea2 = ea->getDerived(PA_GF2());
  
  std::vector<GF2X> slots;
  ea2.decode(slots, encoded);
  
  out.clear();
 
  for (int j = 0; j < 8; ++j) {
#if DEBUG_ENCODE
        cout << "{";
#endif

    Vec4 tmp;
    for (int i = 0; i < 4; ++i)
    {
        int16_t ival = 0;
        BytesFromGF2X((unsigned char*)&ival, slots[(j * 4) + i], 2);
        
#if DEBUG_ENCODE
        cout << ival << ",";
#endif
        tmp[i] = ival;
        
    }
    out.push_back(tmp);
    
#if DEBUG_ENCODE
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
std::vector<Ctxt> ctxt_ylookup;
std::vector<Ctxt> ctxt_xlookup;

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
  
  vertexdata.SetLength(4); // Calculate this
  
  {
    // Idea verticies should be interleaved A0...A7 (slots), B0...B7 (slots) etc.
    EncodeState tmp(ea, vertexdata);
    tmp.EncodeVector({0,0,0,0}); // Bottom Left
    tmp.FinishSlots(); // Interleave all bottom lefts
    tmp.EncodeVector({256,256,0,0}); // Top Middle
    tmp.FinishSlots(); // Interleave all top middle
    tmp.EncodeVector({512,0,0,0}); // Bottom Right
    tmp.FinishSlots(); // Interleave all bottom rights
    tmp.EncodeVector({32767,0,32767,32767}); // Per triangle color -- change to per vertex
    tmp.FinishSlots(); // Interleave all Colors
  
    // encrypt the trangle verticies
    for (int i = 0; i < vertexdata.length(); ++i) {
      ctxt_vertexdata.push_back(Ctxt(*publicKey));
      publicKey->Encrypt(ctxt_vertexdata.back(), vertexdata[i]);
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
             cout << j << " - Fail!!\n";
          }
      }
    }
  }
  
  height = 64;
  width = 64;
  pixelset_size = ea->size() / 4;
  width_pixelset = width / pixelset_size;
  
  // Constants just clear color for now
  // Convert to ctxt_constants
  {
    cout << "Input Constants\n";
    Vec<ZZX> constants;
    constants.SetLength(1);
    EncodeState tmp(ea, constants);
    for (int i = 0; i < pixelset_size; ++i) {
      tmp.EncodeVector({0,0,0,32767});
    }
    ctxt_clearcolor.push_back(Ctxt(*publicKey));
    publicKey->Encrypt(ctxt_clearcolor[0], constants[0]);
  
    for (long i = 0; i < height * width_pixelset; ++i) {
      ctxt_framebuffer.push_back(ctxt_clearcolor[0]);
    }
  }
  
  // Create x/y Lookup Table - find a better way to combine these
  // Do constants need to be encrypted?
  int16_t ps_size = pixelset_size;
  {
    cout << "Y LUT\n";
    Vec<ZZX> ylookup;
    ylookup.SetLength(height);
    
    EncodeState tmp(ea, ylookup);
    for (int16_t y = 0; y < height; ++y) {
      for (int16_t _x = 0; _x < ps_size; ++_x) {
        tmp.EncodeVector({0, y, 0, 0});
      }
    }
    for (int i = 0; i < ylookup.length(); ++i) {
      ctxt_ylookup.push_back(Ctxt(*publicKey));
      publicKey->Encrypt(ctxt_ylookup.back(), ylookup[i]);
    }
  }
  
  {
    cout << "X LUT\n";
    Vec<ZZX> xlookup;
    xlookup.SetLength(width_pixelset);
    
    EncodeState tmp(ea, xlookup);
    for (int16_t x_ = 0; x_ < width_pixelset; ++x_) {
      for (int16_t _x = 0; _x < ps_size; ++_x) {
        int16_t x = (x_ * ps_size) + _x;
        tmp.EncodeVector({x, 0, 0, 0});
      }
    }
    for (int i = 0; i < xlookup.length(); ++i) {
      ctxt_xlookup.push_back(Ctxt(*publicKey));
      publicKey->Encrypt(ctxt_xlookup.back(), xlookup[i]);
    }
  }
    
  // vertex shader rotate matrix
  //transform = shared_ptr<PlaintextMatrixBaseInterface>(buildTransformMatrix(*ea));
}

};

void DuplicateSlots(State& state, const Ctxt& in, Ctxt& out)
{
    const EncryptedArrayDerived<PA_GF2>& ea2 = state.ea->getDerived(PA_GF2());
    out = state.ctxt_clearcolor[0];
    for (int i = 0; i < ea2.size() / 4; ++i) {
        out += in;
        ea2.shift(out, 4);
    }
#if 1
    // test decryption of duplicated slots
    ZZX out1;
    state.secretKey->Decrypt(out1, out);

    vector<Vec4> out2;
    DecodeVectors(state.ea, out1, out2);
    for (int j = 0; j < out2.size(); ++j) {
        if (out2[0] != out2[j]) {
            cout << j << " - Fail!!\n";
        }
    }
#endif
}

void DiscardPoint(State& state, Ctxt& selector,
 const Ctxt& ba, const Ctxt& cb, const Ctxt& ac,
 const Ctxt& a, const Ctxt& b, const Ctxt& c, int x_, int y)
{
    selector = state.ctxt_clearcolor[0]; // Change to 0
    Ctxt p = state.ctxt_ylookup[y];
    p *= state.ctxt_xlookup[x_];
    
#if 1
    // test decryption of the x y co-ord using lookup table
    ZZX p1;
    state.secretKey->Decrypt(p1, p);

    vector<Vec4> p2;
    vector<Vec4> p3;
    DecodeVectors(state.ea, p1, p2);
    for (int16_t _x = 0; _x < state.pixelset_size; ++_x) {
        int16_t x = (x_ * state.pixelset_size) + _x;
        int16_t y1 = y;
        p3.push_back({x, y1, 0, 0});
    }
    assert(p2.size() == p3.size());
    for (int j = 0; j < p2.size(); ++j) {
        if (p2[j] != p3[j]) {
            cout << j << " - Fail!!\n";
        }
    }
#endif
    
    // For now duplicate slot 0 to all slots
    // Later on when there are more than 8 triangles replan what todo here.
    
    // for ba, cb, and ac
    /*DuplicateSlots(state, ba, ba_tmp);
    DuplicateSlots(state, ba, a_tmp);
    tmp -= a_tmp; // p.X - a.x, p.y - a.y, p.z - a.z;
    //swizzle;
    //
    state.ea->rotate(p2, long k);
      
    ba_tmp *= p2; // (b.x - a.x) * (p.y - a.y), (b.y - a.y) * (p.x - a.x)...
    selector *= tmp;*/
}

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
    
    for (int i = 0; i < state.num_triangles; ++i) {
      // Begin Vertex Shader {
      //   mat_mul(*state.ea, state.ctxt_vertexdata[1], *transform);  // transform the triangles

      Ctxt color_out = state.ctxt_vertexdata[(i * 4) + 3];
      // } End Vertex Shader
      
      Ctxt a = state.ctxt_vertexdata[(i * 4) + 0];
      Ctxt b = state.ctxt_vertexdata[(i * 4) + 1];
      Ctxt c = state.ctxt_vertexdata[(i * 4) + 2];
      
      // (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x)
      //      p1            p2            p1            p2
      
      // a-b, b-c, c-a
      
      // Contains paralized data for ea.size() / 4 triangles
      Ctxt ba = b; // b
      ba -= a; // b.x - a.x, b.y - a.y, b.z - a.z;
      Ctxt cb = c; // b
      cb -= b; // c.x - b.x, c.y - b.y, c.z - b.z;
      Ctxt ac = a; // b
      ac -= c; // a.x - c.x, a.y - c.y, a.z - c.z;
      
      // For each pixel fill/rasterize triangle
      for (int y = 0; y < state.height; ++y) {
        for (int x_ = 0; x_ < state.width_pixelset; ++x_) {
          Ctxt discard(*state.publicKey);
          DiscardPoint(state, discard, ba, cb, ac, a, b, c, x_, y);
        
          const int ps = (y * state.width_pixelset) + x_;
          //cout << "PixelSet " << ps << "\n";
          
          Ctxt color_in = color_out;
          color_in *= discard;
          
          // Begin Fragment Shader {
          state.ctxt_framebuffer[ps] += color_in;  // blend current pixel with new
          // } End Fragment Shader
        }
      }
    }
    
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
