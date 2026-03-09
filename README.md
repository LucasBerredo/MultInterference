# MultInterface

Program and accompanying poster to calculate optical interference on multilayer materials, along with some interesting use cases.

Copyright Lucas Berredo, Marcos Vázquez, 2026, licenced under the EUPL.

### Code explanation (flow diagram)

** For individual file explanations, check subdirectories **

```mermaid

graph TD
    subgraph Core Physics Engine
        MI["<b>MultInterference.py</b><br/><i>def MultInterference(...)</i><br/>Calculates optical interference<br/>using transfer matrices"]
    end

    subgraph Application Modules
        BR["<b>Bragg.py</b><br/><i>def bragg(...)</i><br/>Sets up indices/thicknesses<br/>for Bragg Mirrors"]
        FP["<b>Fabry_perot.py</b><br/><i>def fabry_perot(...)</i><br/>Sets up indices/thicknesses<br/>for Fabry-Pérot Cavities"]
    end

    subgraph Execution & Visualization
        EX{"<b>examples.ipynb</b><br/><i>Jupyter Notebook</i><br/>Demonstrates usage and<br/>displays plots"}
    end

    %% Define relationships
    MI -->|Imported by| BR
    MI -->|Imported by| FP
    
    BR -->|Imported & Executed by| EX
    FP -->|Imported & Executed by| EX

    %% Styling
    classDef core fill:#e1f5fe,stroke:#311b92,stroke-width:2px;
    classDef app fill:#f3e5f5,stroke:#4a148c,stroke-width:2px;
    classDef exec fill:#fff3e0,stroke:#880e4f,stroke-width:2px;
    
    class MI core;
    class BR,FP app;
    class EX exec;

```



