//! One module per subcommand. Each exposes `run(args) -> Result<()>`, where `args` are the
//! arguments after the subcommand name.

pub mod intersect;
pub mod merge;
pub mod sort;
pub mod subtract;
pub mod window;
