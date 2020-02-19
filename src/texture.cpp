#include "texture.h"
#include "color.h"

#include <assert.h>
#include <iostream>
#include <algorithm>

using namespace std;

namespace CMU462 {

inline void uint8_to_float( float dst[4], unsigned char* src ) {
  uint8_t* src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8( unsigned char* dst, float src[4] ) {
  uint8_t* dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[0])));
  dst_uint8[1] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[1])));
  dst_uint8[2] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[2])));
  dst_uint8[3] = (uint8_t) ( 255.f * max( 0.0f, min( 1.0f, src[3])));
}

void Sampler2DImp::generate_mips(Texture& tex, int startLevel) {

  // NOTE: 
  // This starter code allocates the mip levels and generates a level 
  // map by filling each level with placeholder data in the form of a 
  // color that differs from its neighbours'. You should instead fill
  // with the correct data!

  // Task 7: Implement this

  // check start level
  if ( startLevel >= tex.mipmap.size() ) {
    std::cerr << "Invalid start level"; 
  }

  // allocate sublevels
  int baseWidth  = tex.mipmap[startLevel].width;
  int baseHeight = tex.mipmap[startLevel].height;
  int numSubLevels = (int)(log2f( (float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  tex.mipmap.resize(startLevel + numSubLevels + 1);

  int width  = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel& level = tex.mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width  = max( 1, width  / 2); assert(width  > 0);
    height = max( 1, height / 2); assert(height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);

  }

  // fill all 0 sub levels with interchanging colors (JUST AS A PLACEHOLDER)
  Color colors[3] = { Color(1,0,0,1), Color(0,1,0,1), Color(0,0,1,1) };
  for(size_t i = 1; i < tex.mipmap.size(); ++i) {

    Color c = colors[i % 3];
    MipLevel& mip = tex.mipmap[i];

    for(size_t i = 0; i < 4 * mip.width * mip.height; i += 4) {
      float_to_uint8( &mip.texels[i], &c.r );
    }
  }

}

Color Sampler2DImp::sample_nearest(Texture& tex, 
                                   float u, float v, 
                                   int level) {

    // Task 6: Implement nearest neighbour interpolation

    // return magenta for invalid level
    if (level < 0 || level > log2f(tex.width) || level > log2f(tex.height))
        return Color(1, 0, 1, 1);
    
    MipLevel& mip = tex.mipmap[level];

    size_t sx = (int)min((size_t)floor(u * (mip.width)), mip.width);
    size_t sy = (int)min((size_t)floor(v * (mip.height)), mip.height);

    float r = mip.texels[4L * (sx + sy * mip.height)] / 255;
    float g = mip.texels[4L * (sx + sy * mip.height) + 1] / 255;
    float b = mip.texels[4L * (sx + sy * mip.height) + 2] / 255;
    float a = mip.texels[4L * (sx + sy * mip.height) + 3] / 255;

    return Color(r, g, b, a);

}

Color Sampler2DImp::sample_bilinear(Texture& tex, 
                                    float u, float v, 
                                    int level) {
  
    // Task 6: Implement bilinear filtering

    // return magenta for invalid level
    if (level < 0 || level > log2f(tex.width) || level > log2f(tex.height))
        return Color(1, 0, 1, 1);
    
    MipLevel& mip = tex.mipmap[level];

    size_t w = mip.width;
    size_t h = mip.height;

    size_t sx = (size_t)floor(u * w);
    size_t sy = (size_t)floor(v * h);

    float x_weight = (u * w - floor(u * w)) / w;
    float y_weight = (u * h - floor(u * h)) / h;    

    Color top_left = sample_nearest(tex, u, v, level);
    Color top_right = sample_nearest(tex, u + 1/w, v, level);
    Color bot_left = sample_nearest(tex, u, v + 1/h, level);
    Color bot_right = sample_nearest(tex, u + 1/w, v + 1/h, level);
    
    float top_r = top_left.r * x_weight + top_right.r * (1 - x_weight);
    float top_g = top_left.g * x_weight + top_right.g * (1 - x_weight);
    float top_b = top_left.b * x_weight + top_right.b * (1 - x_weight);
    float top_a = top_left.a * x_weight + top_right.a * (1 - x_weight);

    float bot_r = bot_left.r * x_weight + bot_right.r * (1 - x_weight);
    float bot_g = bot_left.g * x_weight + bot_right.g * (1 - x_weight);
    float bot_b = bot_left.b * x_weight + bot_right.b * (1 - x_weight);
    float bot_a = bot_left.a * x_weight + bot_right.a * (1 - x_weight);

    return Color(top_r * y_weight + bot_r * (1 - y_weight),
        top_g * y_weight + bot_g * (1 - y_weight), 
        top_b * y_weight + bot_b * (1 - y_weight), 
        top_a * y_weight + bot_a * (1 - y_weight));

}

Color Sampler2DImp::sample_trilinear(Texture& tex, 
                                     float u, float v, 
                                     float u_scale, float v_scale) {

    // Task 7: Implement trilinear filtering

    // return magenta for invalid level
    return Color(1, 0, 1, 1);

}

} // namespace CMU462
