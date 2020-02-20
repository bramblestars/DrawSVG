#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"
#include "texture.h"

using namespace std;

namespace CMU462 {


// Implements SoftwareRenderer //

void SoftwareRendererImp::draw_svg( SVG& svg ) {

  // set top level transformation
  transformation = svg_2_screen;

  // draw all elements
  for ( size_t i = 0; i < svg.elements.size(); ++i ) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  Vector2D a = transform(Vector2D(    0    ,     0    )); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width,     0    )); b.x++; b.y--;
  Vector2D c = transform(Vector2D(    0    ,svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width,svg.height)); d.x++; d.y++;

  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to render target
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;

  w = this->sample_rate * this->target_w;
  h = this->sample_rate * this->target_h;

  sample_buffer.resize(((size_t) 4) * w * h);

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

    // Task 4: 
    // You may want to modify this for supersampling support

    this->render_target = render_target;
    this->target_w = width;
    this->target_h = height;

    w = this->sample_rate * this->target_w;
    h = this->sample_rate * this->target_h;

    sample_buffer.resize(((size_t)4) * w * h);

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

    // Task 5 (part 1):
    // Modify this to implement the transformation stack
      
    Matrix3x3 elem_transform = element->transform;
    Matrix3x3 proj_transform = this->transformation;

    transformation = proj_transform * elem_transform;

    switch(element->type) {
        case POINT:
            draw_point(static_cast<Point&>(*element));
            break;
        case LINE:
            draw_line(static_cast<Line&>(*element));
            break;
        case POLYLINE:
            draw_polyline(static_cast<Polyline&>(*element));
            break;
        case RECT:
            draw_rect(static_cast<Rect&>(*element));
            break;
        case POLYGON:
            draw_polygon(static_cast<Polygon&>(*element));
            break;
        case ELLIPSE:
            draw_ellipse(static_cast<Ellipse&>(*element));
            break;
        case IMAGE:
            draw_image(static_cast<Image&>(*element));
            break;
        case GROUP:
            draw_group(static_cast<Group&>(*element));
            break;
        default:
          break;
      }

    this->transformation = proj_transform;

}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );

}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Extra credit 

}

void SoftwareRendererImp::draw_image( Image& image ) {

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

    int rate = this->sample_rate;

    // fill in the nearest pixel
    int sx = (int) floor(x);
    int sy = (int) floor(y);

    sx *= rate;
    sy *= rate;

    // check bounds
    if ( sx < 0 || sx >= w ) return;
    if ( sy < 0 || sy >= h ) return;

    // fill sample - NOT doing alpha blending!
    fill_sample(sx, sy, color);

}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {

    int sx0 = (int)floor(x0);
    int sy0 = (int)floor(y0);
    int sx1 = (int)floor(x1);
    int sy1 = (int)floor(y1);

    //if vertical line
    if (sx0 == sx1) {

        for (int i = min(sy0, sy1); i <= max(sy0, sy1); i++) {
            rasterize_point(sx0, i, color);
        }

        return;

    }
   
    float slope = (y1 - y0) / (x1 - x0);
    int x_start, y_start;
    float x, y;

    if (-1 <= slope && slope <= 1) {

         if (sx0 < sx1) {
             x_start = sx0;
             y_start = sy0;
         }

         else { //sx0 >= sx1
             x_start = sx1;
             y_start = sy1;
         }

         y = y_start == sy0? y0 : y1;

         for (int i = x_start; i <= max(sx0, sx1); i++) {
             rasterize_point(i, y, color);
             y += slope;
         }
    }

    else {
         if (sy0 < sy1) {
             x_start = sx0;
             y_start = sy0;
         }

         else { //sy0 >= sy1
             x_start = sx1;
             y_start = sy1;
         }

         x = x_start == sx0? x0 : x1;
         
         for (int i = y_start; i <= max(sy0, sy1); i++) {
             rasterize_point(x, i, color);
             x += (1 / slope);
         }

    }

}

float signed_area(float x0, float y0, float x1, float y1, float x2, float y2) {
    return (x0 * y1 - x1 * y0) - (x0 * y2 - x2 * y0) + (x1 * y2 - x2 * y1);
}


void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {

    float rate = (float)this->sample_rate;

    float bottom = min(y0, min(y1, y2)) * rate;
    float top = max(y0, max(y1, y2)) * rate;
    float left = min(x0, min(x1, x2)) * rate;
    float right = max(x0, max(x1, x2)) * rate;

    float a1, a2, a3;
    float x_test, y_test;

    for (float x = floor(left); x <= right; x++) {
        for (float y = floor(bottom); y <= top; y++) {

            x_test = x + 0.5;
            y_test = y + 0.5;

            a1 = signed_area(x_test, y_test, x0, y0, x1, y1);
            a2 = signed_area(x_test, y_test, x1, y1, x2, y2);
            a3 = signed_area(x_test, y_test, x2, y2, x0, y0);

            if ((a1 >= 0 && a2 >= 0 && a3 >= 0) || (a1 <= 0 && a2 <= 0 && a3 <= 0)) {
                rasterize_point(x/rate, y/rate, color);
            }

        }
    }

}

void SoftwareRendererImp::rasterize_image(float x0, float y0,
    float x1, float y1,
    Texture& tex) {
    // Task 6: 
    // Implement image rasterization

    int rate = this->sample_rate;

    size_t sx0 = (size_t)floor(x0 * rate);
    size_t sx1 = (size_t)floor(x1 * rate);
    size_t sy0 = (size_t)floor(y0 * rate);
    size_t sy1 = (size_t)floor(y1 * rate);

    float u, v, du, dv;

    for (size_t x = sx0; x < sx1; x++) {
        for (size_t y = sy0; y < sy1; y++) {
            u = (float)(x - sx0) / (float)(sx1 - sx0);
            v = (float)(y - sy0) / (float)(sy1 - sy0);
            
            du = 1 / (sx1 - sx0);
            dv = 1 / (sx1 - sx0);

            fill_sample(x, y, sampler->sample_trilinear(tex, u, v, du, dv));

        }
    }

}

void SoftwareRendererImp::fill_sample(int sx, int sy, const Color& c) {
    sample_buffer[4L * ((size_t) sx + (size_t) sy * w)    ] = (uint8_t) (c.r * 255);
    sample_buffer[4L * ((size_t) sx + (size_t) sy * w) + 1] = (uint8_t) (c.g * 255);
    sample_buffer[4L * ((size_t) sx + (size_t) sy * w) + 2] = (uint8_t) (c.b * 255);
    sample_buffer[4L * ((size_t) sx + (size_t) sy * w) + 3] = (uint8_t) (c.a * 255);
}

void SoftwareRendererImp::fill_pixel(int x, int y, const Color& c) {
    this->render_target[4 * (x + y * this->target_w)    ] = (uint8_t) (c.r * 255);
    this->render_target[4 * (x + y * this->target_w) + 1] = (uint8_t) (c.g * 255);
    this->render_target[4 * (x + y * this->target_w) + 2] = (uint8_t) (c.b * 255);
    this->render_target[4 * (x + y * this->target_w) + 3] = (uint8_t) (c.a * 255);
}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".

    float r, b, g, a;
    size_t rate = this->sample_rate;

    for (size_t x = 0; x < target_w; x++) {
        for (size_t y = 0; y < target_h; y++) {
            
            r = 0;
            b = 0;
            g = 0;
            a = 0;

            for (size_t sx = x * rate; sx < x * rate + rate; sx++) {
                for (size_t sy = y * rate; sy < y * rate + rate; sy++) {
                    r += (float) sample_buffer[4L * (sx + sy * w)];
                    g += (float) sample_buffer[4L * (sx + sy * w) + 1L];
                    b += (float) sample_buffer[4L * (sx + sy * w) + 2L];
                    a += (float) sample_buffer[4L * (sx + sy * w) + 3L];

                    sample_buffer[4L * (sx + sy * w)] = 255;
                    sample_buffer[4L * (sx + sy * w) + 1L] = 255;
                    sample_buffer[4L * (sx + sy * w) + 2L] = 255;
                    sample_buffer[4L * (sx + sy * w) + 3L] = 255;
                }
            }

            r = r / (rate * rate) / 255;
            g = g / (rate * rate) / 255;
            b = b / (rate * rate) / 255;
            a = a / (rate * rate) / 255;

            fill_pixel(x, y, Color(r, g, b, a));

        }
    }


}


} // namespace CMU462
