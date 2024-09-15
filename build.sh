HOST=$(rustc -vV | rg host | sed 's/^host: //g')
RUSTFLAGS=-Zlocation-detail=none \
  cargo +nightly build \
  -Z build-std=std,panic_abort \
  -Z build-std-features="optimize_for_size" \
  --target $HOST \
  --release
