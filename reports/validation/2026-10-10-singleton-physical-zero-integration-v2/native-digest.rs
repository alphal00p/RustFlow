//! Validation-only adapter to the existing native BLAKE3 implementation.
use std::io::Read;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args_os();
    args.next();
    let path = args.next().ok_or("one file path is required")?;
    if args.next().is_some() {
        return Err("one file path is required".into());
    }
    let mut input = std::fs::File::open(path)?;
    let mut hasher = blake3::Hasher::new();
    let mut buffer = [0_u8; 65536];
    let mut bytes: u64 = 0;
    loop {
        let count = input.read(&mut buffer)?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
        bytes = bytes.checked_add(count as u64).ok_or("file too large")?;
    }
    println!("{} {bytes}", hasher.finalize().to_hex());
    Ok(())
}
