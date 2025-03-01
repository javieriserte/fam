use std::{env, path::Path};
use serde::Deserialize;
use serde;
use home::home_dir;


#[derive(Clone, Debug, Deserialize)]
pub struct Config {
    #[serde(default)]
    pub locale: Locale
}

fn en() -> String { "en".to_string() }

#[derive(Clone, Debug, Deserialize)]
pub struct Locale {
    #[serde(default = "en")]
    pub lang: String
}

impl Default for Locale {
    fn default() -> Self {
        Locale {
            lang: "en".to_string()
        }
    }
}

pub fn get_locale_from_env() -> String {
    let locale = env::var("LANG")
        .or_else(|_| env::var("LC_ALL"))
        .unwrap_or_else(|_| "en-US".to_string());
    locale
}

pub fn default_config() -> Config {
    Config {
        locale: Locale{
            lang: get_locale_from_env()
        }
    }
}

pub fn get_config_from_file() -> Config {
    let config = default_config();
    let home = home_dir()
        .map(
            |x| x
                .join(".config")
                .join("fam")
                .join("config.toml")
        );
    match home {
        Some(path) => {
            if path.exists() {
                let config_file = path.to_str().unwrap();
                let config_path = Path::new(config_file);
                std::fs::read_to_string(config_path)
                    .map_err(
                        |_|"[WARN] Could not read config file, using default config"
                    )
                    .and_then(
                        |x|
                            toml::from_str::<Config>(&x)
                            .map_err(|_| "Could not parse config file")
                    )
                    .map(
                        |x| {
                            let mut c = config.clone();
                            c.locale.lang = x.locale.lang;
                            c
                        }
                    )
                    .unwrap_or_else(
                        |x| {
                            eprintln!("{}", x);
                            config
                        }
                    )
            } else {
                config
            }
        },
        None => config
    }
}

pub fn get_config() -> Config {
    get_config_from_file()
}
