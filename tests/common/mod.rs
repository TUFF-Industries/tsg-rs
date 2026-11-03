//! Builds small synthetic TSG pairs in a temporary directory. No real data.

#![allow(dead_code)]

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

/// TSG's scalar null marker.
pub const NULL: f32 = f32::from_bits(0xFF7F_FFFE);

/// A synthetic NIR/SWIR pair: 3 samples x 5 bands on a 1300-2500 nm grid.
/// Scalar columns: 0 `Depth (m)` (plain), 1 `Grp1 sTSAS` (class 1),
/// 2 `Min1 sTSAS` (class 2), 3 `Grp1 sTSAV` (class 0).
pub struct Synth {
    pub header: String,
    pub spectra: Vec<f32>,
    pub scalars: Vec<f32>,
}

pub const SAMPLES: usize = 3;
pub const BANDS: usize = 5;

impl Synth {
    pub fn new() -> Self {
        let header = "\
[description]
SYNTHETIC TEST FILE

[structure]
datatype = float
dataheader = 0

[coordinates]
firstsample = 1
lastsample = 3
firstline = 1
lastline = 1
firstband = 1
lastband = 5

[reflibs]
system = SYS_TEST_SWIR

[band headers]
0:Depth (m);10.000000;10.200000;2;-1;0
1:Grp1 sTSAS;1.000000;3.000000;2;1;1;1;0
2:Min1 sTSAS;1.000000;3.000000;2;2;1;0;0
3:Grp1 sTSAV;1.000000;3.000000;2;0;2;1;0

[sample headers]
0:SYN_0001_1  T=0001 L=1 P=1 D=10.000000 X=1.0 H=SYN
1:SYN_0001_2  T=0001 L=1 P=2 D=10.100000 X=2.0 H=SYN
2:SYN_0001_3  T=0001 L=1 P=3 D=10.200000 X=3.0 H=SYN

[class 0]
name = S_VNIR_TEST Groups
max = 2
0:OXIDE
1:INVALID

[class 1]
name = S_SWIR_TSA705 Groups
max = 3
colours = 1 2 3
0:KAOLIN
1:WHITE-MICA
2:INVALID

[class 2]
name = S_SWIR_TSA705 Minerals
max = 2
0:Kaolinite-WX
1:Muscovite

[class 3]
name = Domain
max = 1
0:Default

[wavelength specs]
1300.000000 2500.000000 nm
"
        .to_string();
        let spectra = (0..SAMPLES * BANDS).map(|i| 0.1 + i as f32 / 7.0).collect();
        #[rustfmt::skip]
        let scalars = vec![
            10.0, 0.0, 1.0, 0.0, NULL,
            10.1, 1.0, 1.0, 1.0, NULL,
            10.2, NULL, NULL, NULL, 0.5,
        ];
        Self {
            header,
            spectra,
            scalars,
        }
    }

    pub fn bip(&self) -> Vec<u8> {
        self.spectra
            .iter()
            .chain(&self.scalars)
            .flat_map(|v| v.to_le_bytes())
            .collect()
    }

    /// Write `<dir>/SYN_tsg.tsg` and `.bip`; returns the `.tsg` path.
    pub fn write(&self, dir: &Path) -> PathBuf {
        let tsg = dir.join("SYN_tsg.tsg");
        std::fs::write(&tsg, &self.header).unwrap();
        std::fs::write(dir.join("SYN_tsg.bip"), self.bip()).unwrap();
        tsg
    }
}

/// A fresh empty directory under the system temp dir.
pub fn temp_dir(tag: &str) -> PathBuf {
    static N: AtomicUsize = AtomicUsize::new(0);
    let d = std::env::temp_dir().join(format!(
        "tsg-rs-test-{}-{tag}-{}",
        std::process::id(),
        N.fetch_add(1, Ordering::SeqCst)
    ));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}
