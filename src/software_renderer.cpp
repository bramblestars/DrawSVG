#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

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

}

void SoftwareRendererImp::set_render_target( unsigned char* render_target,
                                             size_t width, size_t height ) {

  // Task 4: 
  // You may want to modify this for supersampling support
  this->render_target = render_target;
  this->target_w = width;
  this->target_h = height;

}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

  // Task 5 (part 1):
  // Modify this to implement the transformation stack

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

  // fill in the nearest pixel
  int sx = (int) floor(x);
  int sy = (int) floor(y);

  // check bounds
  if ( sx < 0 || sx >= target_w ) return;
  if ( sy < 0 || sy >= target_h ) return;

  // fill sample - NOT doing alpha blending!
  render_target[4 * (sx + sy * target_w)    ] = (uint8_t) (color.r * 255);
  render_target[4 * (sx + sy * target_w) + 1] = (uint8_t) (color.g * 255);
  render_target[4 * (sx + sy * target_w) + 2] = (uint8_t) (color.b * 255);
  render_target[4 * (sx + sy * target_w) + 3] = (uint8_t) (color.a * 255);

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

float cross_product(float x0, float y0,
                   float x1, float y1) {
    return (x0 * y1) - (x1 * y0);
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
    /*
    int sx0 = (int)floor(x0); int sy0 = (int)floor(y0);
    int sx1 = (int)floor(x1); int sy1 = (int)floor(y1);
    int sx2 = (int)floor(x2); int sy2 = (int)floor(y2);

    int box_left = min(sx0, max(sx1, sx2));
    int box_right = max(sx0, max(sx1, sx2));
    int box_bottom = min(sy0, max(sy1, sy2));
    int box_top = max(sy0, max(sy1, sy2));

    for (int x = box_left; x <= box_right; x++) {
        for (int y = box_bottom; y <= box_top; y++) {

            float det = cross_product(x1 - x0, y1 - y0, x2 - x0, y2 - y0);
            float s = cross_product(x - x0, y - y0, x2 - x0, y2 - y0) / det;
            float t1 = cross_product(x - x0, y - y0, x1 - x0, y1 - y0) / det;
            float t2 = cross_product(x1 - x0, y1 - y0, x - x0, y - y0) / det;

            if ((s >= 0 && t1 >= 0 && s + t1 <= 1) || 
                (s >= 0 && t2 >= 0 && s + t2 <= 1)) {
                rasterize_point(x, y, color);
            }

        }
    }
    */

    float slope1, slope2;
    float curr_left, curr_right;

    // Sort vertices from top to bottom 
    float bottom_x, bottom_y, middle_x, middle_y, top_x, top_y;

    bottom_y = min(y0, min(y1, y2));
    middle_y = min(min(max(y0, y1), max(y1, y2)), max(y0, y2));
    top_y = max(y0, max(y1, y2));

    if (bottom_y == y0) {
        bottom_x = x0;

        if (middle_y == y1) {
            middle_x = x1;
            top_x = x2;
        }
        else {
            middle_x = x2;
            top_x = x1;
        }
    }

    else if (bottom_y == y1) {
        bottom_x = x1;

        if (middle_y == y0) {
            middle_x = x0;
            top_x = x2;
        }
        else {
            middle_x = x2;
            top_x = x0;
        }
    }

    else {
        bottom_x = x2;

        if (middle_y == y0) {
            middle_x = x0;
            top_x = x1;
        }
        else {
            middle_x = x1;
            top_x = x0;
        }
    }

    // if bottom edge is flat
    if (floor(bottom_y) == floor(middle_y)) {

        curr_left = min(bottom_x, middle_x);
        curr_right = max(bottom_x, middle_x);

        slope1 = (top_x - curr_left) / (top_y - bottom_y);
        slope2 = (top_x - curr_right) / (top_y - bottom_y);

        //draw horizontal lines starting at the bottom and going to the top
        for (int y = (int) floor(bottom_y); y <= (int) floor(top_y); y++) {
            rasterize_line(curr_left, y, curr_right, y, color);
            curr_left += slope2;
            curr_right += slope1;
        }

    }

    // if top edge is flat
    else if (floor(middle_y) == floor(top_y)) {

        curr_left = min(middle_x, top_x);
        curr_right = max(middle_x, top_x);

        slope1 = (curr_left - bottom_x) / (top_y - bottom_y);
        slope2 = (curr_right - bottom_x) / (top_y - bottom_y);

        //draw horizontal lines starting at the bottom and going to the top
        for (int y = (int) floor(top_y); y >= (int)floor(bottom_y); y--) {
            rasterize_line(curr_left, y, curr_right, y, color);
            curr_left += slope2;
            curr_right += slope1;
        }
    }

    // otherwise, split triangle into a bottom-flat triangle and a top-flat triangle
    
    else {
        float x = bottom_x + (middle_y - bottom_y) / (top_y - bottom_y) * (top_x - bottom_x);
        float y = middle_y;

        rasterize_triangle(middle_x, middle_y, x, y, top_x, top_y, color);
        rasterize_triangle(middle_x, middle_y, x, y, bottom_x, bottom_y, color);
    }
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 6: 
  // Implement image rasterization

}

// resolve samples to render target
void SoftwareRendererImp::resolve( void ) {

  // Task 4: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 4".
  return;

}


} // namespace CMU462
