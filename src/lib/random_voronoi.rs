use std::{path::Path};
use rand::{Rng, prelude::SliceRandom, thread_rng};
use graphics_buffer::RenderBuffer;


pub fn random_color() -> [f32; 4] {
    [
        rand::thread_rng().gen_range(0..=100) as f32/100f32,
        rand::thread_rng().gen_range(0..=100) as f32/100f32,
        rand::thread_rng().gen_range(0..=100) as f32/100f32,
        1.0
    ]
}

pub fn random_point() -> [u32; 2] {
    [
        rand::thread_rng().gen_range(0..=1000),
        rand::thread_rng().gen_range(0..=1000)
    ]
}

pub fn build_points(npoints:usize) -> Vec<([u32;2],[f32;4])> {
    (0..npoints)
        .map(|_| (random_point(), random_color()))
        .collect()
}

pub fn dist(p1: &[u32;2], p2: &[u32;2]) -> f64 {
    (((p1[0] as i32-p2[0] as i32) as f64).powf(2f64)+((p1[1] as i32-p2[1] as i32) as f64).powf(2f64)).sqrt()
}

pub fn select_close_point(p: &[u32;2], points: &Vec<([u32;2],[f32;4])>) -> [f32;4]{
    let many = 100;
    let mut rng = thread_rng();
    let mut shuffled = (0..points.len()).into_iter().collect::<Vec<usize>>();
    shuffled.shuffle(&mut rng);
    shuffled = shuffled.into_iter().take(many).collect();
    let mut closest: Option<usize> = None;
    let mut closest_dist = 1000f64;
    shuffled.iter().for_each(
        |x| {
            match closest {
                None => {
                    closest = Some(*x);
                    closest_dist = dist(p, &points[*x].0);
                },
                Some(_) => {
                    let c_dist =  dist(p, &points[*x].0);
                    if c_dist < closest_dist {
                        closest_dist = c_dist;
                        closest = Some(*x)
                    }
                }
            }
        }
    );
    match closest {
        Some(x) => {
            points[x].1
        }
        None => {
            [0.0,0.0,0.0,1.0]
        }
    }
}

pub fn random_voro(outfile: &Path) {
    let width:u32 = 1000;
    let height:u32 = 1000;
    let points = build_points(300);
    let mut buffer = RenderBuffer::new(width, height);

    buffer.clear([1.0, 1.0, 1.0, 1.0]);
    for x in 0..width {
        for y in 0..height {
            let color = select_close_point(&[x, y], &points);
            buffer.set_pixel(x, y , color);
        }
        if x % 100 == 0 {
            println!("rows: {}", x)
        }
    }

    // Save the buffer
    let x = outfile.display().to_string();
    buffer.save(&x[..]).unwrap();
}
#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_random_voro() {
        let outfile = "a.png";
        random_voro(Path::new(outfile))
    }
}