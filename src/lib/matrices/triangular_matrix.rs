use std::ops::{Add, AddAssign};
use std::{error::Error, io::ErrorKind};

/// Declares the `One` trait for `usize`.
///
/// This trait provides a method `one()` that returns the value `1` for
/// different numeric types.
///
/// # Examples
/// ```
/// use fam::matrices::triangular_matrix::One;
/// let one = usize::one();
/// assert_eq!(one, 1);
/// ```
pub trait One { fn one() -> Self ;}
impl One for usize { fn one() -> Self { 1 as usize } }
impl One for i32 { fn one() -> Self { 1 as i32 } }
impl One for u32 { fn one() -> Self { 1 as u32 } }
impl One for f32 { fn one() -> Self { 1 as f32 } }
impl One for f64 { fn one() -> Self { 1 as f64 } }


/// Declares the `Num` trait that represents a numeric type that can be used
/// in the `TriangularMatrix` struct.
pub trait Num:
    Add<Output = Self> +
    AddAssign<Self> +
    One +
    Copy +
    Default +
    Clone {}

impl<T> Num for T
where T:
    Add<Output = T> +
    AddAssign<T> +
    One +
    Copy +
    Default +
    Clone {}


/// Represents a triangular matrix.
///
/// A triangular matrix is a symmetric matrix where all elements below and above
/// the diagonal are identical.
pub struct TriangularMatrix<T> where T: Num
{
    pub data: Vec<T>,
    pub size: usize
}

impl <T> TriangularMatrix<T> where T: Num
{
    /// Creates a new `TriangularMatrix` with a given width.
    pub fn new(width: usize) -> Self {
        let size: usize = width*(width+1)/2;
        Self{
            data: vec![T::default(); size],
            size
        }
    }
    /// Gets the value at a given position.
    pub fn get(&self, x: usize, y: usize) -> Result<T, Box<dyn Error>> {
        let pos = self.xy_to_linear_coord(x, y)?;
        Ok(self.data[pos])
    }
    /// Sets the value at a given position.
    pub fn set(
        &mut self,
        x: usize,
        y: usize,
        value: T
    ) -> Result<(), Box<dyn Error>> {
        let pos = self.xy_to_linear_coord(x, y)?;
        self.data[pos] = value;
        Ok(())
    }
    /// Returns the size of the matrix.
    pub fn size(&self) -> usize {
        self.size
    }
    /// Converts a linear index to a pair of x and y coordinates.
    /// The linear index is the index of the element in the data vector.
    /// Although this method is public, it is not intended to be used outside
    /// of the struct.
    pub fn linear_to_xy_coord(
        &self,
        pos: usize
    ) -> Result<(usize, usize), Box<dyn Error>> {
        if pos < self.size() {
            let v = (8f64*pos as f64 + 1f64).sqrt().floor() as usize;
            let y = (v-1) / 2;
            let x = pos - y*(y+1)/2;
            Ok((x, y))
        } else {
            Err(
                Box::new(
                    std::io::Error::new(
                        ErrorKind::Other,
                        format!("Index out of range")
                    )
                )
            )
        }
    }
    /// Converts a pair of x and y coordinates to a linear index.
    /// The linear index is the index of the element in the data vector.
    /// Although this method is public, it is not intended to be used outside
    /// of the struct.
    pub fn xy_to_linear_coord(
        &self,
        x: usize,
        y: usize
    ) -> Result<usize, Box<dyn Error>> {
        let pos = match y>=x {
            true => y*(y+1)/2+x,
            false => x*(x+1)/2+y
        };
        if pos < self.size() {
            Ok(pos)
        } else {
            Err(
                Box::new(
                    std::io::Error::new(
                        ErrorKind::Other,
                        format!("Index out of range")
                    )
                )
            )
        }
    }
    /// Increments the value at a given position by one.
    pub fn increment(&mut self, x:usize, y:usize) -> Result<(), Box<dyn Error>> {
        let pos = self.xy_to_linear_coord(x, y)?;
        let one = T::one();
        self.data[pos] += one;
        Ok(())
    }
}


#[cfg(test)]
mod test {
    use crate::matrices::TriangularMatrix;
    #[test]
    fn test_new() {
        let tm = TriangularMatrix::<usize>::new(3);
        assert_eq!(tm.size(), 6);
        assert!(tm.data.iter().all(|&x| x == 0));
    }
    #[test]
    fn test_get_set() {
        let mut tm = TriangularMatrix::<usize>::new(3);
        tm.set(0, 1, 1).unwrap();
        assert_eq!(tm.get(0, 1).unwrap(), 1);
        tm.set(1, 0, 2).unwrap();
        assert_eq!(tm.get(0, 1).unwrap(), 2);
        assert_eq!(tm.get(1, 0).unwrap(), 2);
        assert_eq!(tm.get(4, 4).is_err(), true);
    }
    #[test]
    fn test_linear_to_xy_coord() {
        let tm = TriangularMatrix::<usize>::new(3);
        assert_eq!(tm.linear_to_xy_coord(0).unwrap(), (0, 0));
        assert_eq!(tm.linear_to_xy_coord(1).unwrap(), (0, 1));
        assert_eq!(tm.linear_to_xy_coord(2).unwrap(), (1, 1));
        assert_eq!(tm.linear_to_xy_coord(3).unwrap(), (0, 2));
        assert_eq!(tm.linear_to_xy_coord(4).unwrap(), (1, 2));
        assert_eq!(tm.linear_to_xy_coord(5).unwrap(), (2, 2));
        assert_eq!(tm.linear_to_xy_coord(6).is_err(), true);
    }
    #[test]
    fn test_xy_to_linear_coord() {
        let tm = TriangularMatrix::<usize>::new(3);
        assert_eq!(tm.xy_to_linear_coord(0, 0).unwrap(), 0);
        assert_eq!(tm.xy_to_linear_coord(0, 1).unwrap(), 1);
        assert_eq!(tm.xy_to_linear_coord(1, 0).unwrap(), 1);
        assert_eq!(tm.xy_to_linear_coord(1, 1).unwrap(), 2);
        assert_eq!(tm.xy_to_linear_coord(0, 2).unwrap(), 3);
        assert_eq!(tm.xy_to_linear_coord(2, 0).unwrap(), 3);
        assert_eq!(tm.xy_to_linear_coord(1, 2).unwrap(), 4);
        assert_eq!(tm.xy_to_linear_coord(2, 1).unwrap(), 4);
        assert_eq!(tm.xy_to_linear_coord(2, 2).unwrap(), 5);
        assert_eq!(tm.xy_to_linear_coord(0, 3).is_err(), true);
    }
    #[test]
    fn test_increment() {
        let mut tm = TriangularMatrix::<usize>::new(3);
        tm.increment(0, 1).unwrap();
        assert_eq!(tm.get(0, 1).unwrap(), 1);
        tm.increment(0, 1).unwrap();
        assert_eq!(tm.get(0, 1).unwrap(), 2);
        tm.increment(1, 0).unwrap();
        assert_eq!(tm.get(0, 1).unwrap(), 3);
        assert_eq!(tm.get(1, 0).unwrap(), 3);
    }
}