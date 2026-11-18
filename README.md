# John Kaluiji

Personal website: résumé, projects, and engineering simulations that run in the browser.

Plain HTML, CSS and JavaScript with no build step. Each simulation has its own page and a matching `math-*.html` page that works through its equations.

| Page | What it is |
|---|---|
| `index.html` | Landing page |
| `resume.html`, `resume-print.html` | Full résumé and the one-page version |
| `simulations.html` | Index of the simulations |
| `airflow.html` + `wake-model.js` | Navier-Stokes flow past a cylinder (vortex shedding) and ideal flow |
| `convection.html` + `convection-model.js` | Forced convection over a heated aluminum plate |
| `beam.html` + `beam-model.js` | Cantilever beam, Euler-Bernoulli finite elements |
| `servo-arm.html` | Four-joint arm driven by 555 timers, with schematics and PCB |
| `solar-sim.html`, `solar-tracking.html` + `solar-model.js` | Solar tracking simulation and case study |

Run locally: open `index.html`, or serve the folder with `python -m http.server`.
