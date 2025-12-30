//! Circular array operations for genomic data.
//!
//! Handles circular genomes where reads can wrap around from the end
//! back to the beginning (e.g., position 9900-100 in a 10kb genome).

use std::ops::AddAssign;

// ============================================================================
// CIRCULAR ARRAY FUNCTIONS
// ============================================================================

/// Increment values in a circular range [start, end), handling wrap-around.
#[inline]
pub fn increment_circular<T>(arr: &mut [T], start: usize, end: usize, delta: T)
where
    T: AddAssign + Copy,
{
    let len = arr.len();
    let start_mod = start % len;
    let end_mod = end % len;

    if start_mod < end_mod {
        // Fast path: no wrapping
        for pos in start_mod..end_mod {
            arr[pos] += delta;
        }
    } else {
        // Wrapped or edge case
        let read_len = end - start;
        if read_len == 0 {
            return;
        }
        if read_len < len {
            // Simple wrap
            for pos in start_mod..len {
                arr[pos] += delta;
            }
            for pos in 0..end_mod {
                arr[pos] += delta;
            }
        } else {
            // Read >= genome length (rare)
            let full_laps = read_len / len;
            let remaining = read_len % len;

            for _ in 0..full_laps {
                for pos in 0..len {
                    arr[pos] += delta;
                }
            }

            if remaining > 0 {
                let rem_end = (start_mod + remaining) % len;
                if start_mod < rem_end {
                    for pos in start_mod..rem_end {
                        arr[pos] += delta;
                    }
                } else {
                    for pos in start_mod..len {
                        arr[pos] += delta;
                    }
                    for pos in 0..rem_end {
                        arr[pos] += delta;
                    }
                }
            }
        }
    }
}

/// Alias for increment_circular (for compatibility)
#[inline]
pub fn increment_circular_long<T>(arr: &mut [T], start: usize, end: usize, delta: T)
where
    T: AddAssign + Copy,
{
    increment_circular(arr, start, end, delta);
}

/// Increment values in a non-circular range [start, end).
#[inline]
pub fn increment_range<T>(arr: &mut [T], start: usize, end: usize, delta: T)
where
    T: AddAssign + Copy,
{
    for pos in start..end {
        arr[pos] += delta;
    }
}
