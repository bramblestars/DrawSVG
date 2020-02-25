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

  
  Color c;
  float r, g, b, a;
  size_t w, h;

  for (size_t i = 1; i < tex.mipmap.size(); ++i) {

      MipLevel& prev_mip = tex.mipmap[i - 1];
      MipLevel& mip = tex.mipmap[i];
      w = mip.width; h = mip.height;

      for (size_t i = 0; i < 4 * w * h; i += 4) {
          r = 0; g = 0; b = 0; a = 0;

          r += (float)prev_mip.texels[i * 2];
          g += (float)prev_mip.texels[i * 2 + 1];
          b += (float)prev_mip.texels[i * 2 + 2];
          a += (float)prev_mip.texels[i * 2 + 3];

          r += (float)prev_mip.texels[i * 2 + 5];
          g += (float)prev_mip.texels[i * 2 + 6];
          b += (float)prev_mip.texels[i * 2 + 7];
          a += (float)prev_mip.texels[i * 2 + 8];

          r += (float)prev_mip.texels[(i + w * 4) * 2];
          g += (float)prev_mip.texels[(i + w * 4) * 2 + 1];
          b += (float)prev_mip.texels[(i + w * 4) * 2 + 2];
          a += (float)prev_mip.texels[(i + w * 4) * 2 + 3];

          r += (float)prev_mip.texels[(i + w * 4) * 2 + 5];
          g += (float)prev_mip.texels[(i + w * 4) * 2 + 6];
          b += (float)prev_mip.texels[(i + w * 4) * 2 + 7];
          a += (float)prev_mip.texels[(i + w * 4) * 2 + 8];

          r /= 4;
          g /= 4;
          b /= 4;
          a /= 4;

          mip.texels[i] = r;
          mip.texels[i + 1] = g;
          mip.texels[i + 2] = b;
          mip.texels[i + 3] = a;
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

    size_t sx = (size_t)floor(u * (mip.width));
    size_t sy = (size_t)floor(v * (mip.height));

    Color color;

    uint8_to_float(&color.r, &mip.texels[4L * (sx + sy * mip.width)]);

    return color;

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

    float x_weight = 1 - (u * w - floor(u * w)) / (float) w;
    float y_weight = 1 - (u * h - floor(u * h)) / (float) h;    

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
    
    size_t w = tex.width;
    size_t h = tex.height;

    int max_level = min(log2f(w), log2f(h));

    float dx = 1 / w;
    float dy = 1 / h;

    float L = max(u_scale / dx, v_scale / dy);
    float mip_level = log2f(L);
    int level_below = (int)floor(mip_level);
    int level_above = (int)ceil(mip_level);
    
    // return magenta for invalid level
    if (level_below > max_level) {
        return Color(1, 0, 1, 1);
    }

    if (level_below < 0) {
        return sample_bilinear(tex, u, v, 0);
    }

    if (level_above > log2f(w) || level_above > log2f(h)) {
        return sample_bilinear(tex, u, v, max_level);
    }
    
    float level_weight = mip_level - level_below;

    Color above = sample_bilinear(tex, u, v, level_above);
    Color below = sample_bilinear(tex, u, v, level_below);

    float r = level_weight * above.r + (1 - level_weight) * below.r;
    float g = level_weight * above.g + (1 - level_weight) * below.g;
    float b = level_weight * above.b + (1 - level_weight) * below.b;
    float a = level_weight * above.a + (1 - level_weight) * below.a;

    return Color(r, b, g, a);

}

} // namespace CMU462
