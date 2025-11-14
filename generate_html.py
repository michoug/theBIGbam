"""
Simple Flask app that generates an interactive HTML page.

- Connects to an sqlite DB (passed as --db argument).
- The page contains a collapsible control panel with:
  - Sample select (populated from `Sample` table)
  - Contig select (populated from `Contig` table)
  - Checkboxes for each distinct Module found in `Variable.Module`
  - Apply button
- When Apply is pressed the page sends an AJAX request to the server which:
  - looks up variables for the selected modules (from `Variable` table)
  - finds Feature_<variable> tables and queries rows for the chosen Sample and Contig
  - returns a small extract of the data which is displayed in the main panel

This file intentionally keeps dependencies minimal: Flask is required to run the server.

Usage:
    python generate_html.py --db results.db --host 127.0.0.1 --port 5000

Then open http://127.0.0.1:5000/ in your browser.

"""

import argparse
import sqlite3
import json
import html
import os
from flask import Flask, request, jsonify, render_template_string

app = Flask(__name__)
DB_PATH = None

# Simple HTML template using Bootstrap 4 from CDN and a bit of JS
HTML_TEMPLATE = """
<!doctype html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
    <title>Feature viewer (DB-driven)</title>
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.5.2/css/bootstrap.min.css">
    <style>
      body { padding: 1rem; }
      #sidebar { max-width: 320px; }
      #main-panel { margin-left: 0; }
      pre.json { background:#f8f9fa; padding: 1rem; border-radius: 4px; }
    </style>
  </head>
  <body>
    <div class="container-fluid">
      <div class="row">
        <div class="col-md-3" id="sidebar">
          <div class="card">
            <div class="card-header">
              Controls
              <button class="btn btn-sm btn-outline-secondary float-right" id="toggle-btn">Hide</button>
            </div>
            <div class="card-body" id="controls">
              <form id="filter-form">
                <div class="form-group">
                  <label for="sample-select">Sample</label>
                  <select class="form-control" id="sample-select" name="sample">
                    {% for s in samples %}
                      <option value="{{s}}">{{s}}</option>
                    {% endfor %}
                  </select>
                </div>
                <div class="form-group">
                  <label for="contig-select">Contig</label>
                  <select class="form-control" id="contig-select" name="contig">
                    {% for c in contigs %}
                      <option value="{{c}}">{{c}}</option>
                    {% endfor %}
                  </select>
                </div>
                <div class="form-group">
                  <label>Modules</label>
                  <div>
                    {% for m in modules %}
                      <div class="form-check">
                        <input class="form-check-input module-checkbox" type="checkbox" value="{{m}}" id="mod_{{loop.index}}">
                        <label class="form-check-label" for="mod_{{loop.index}}">{{m}}</label>
                      </div>
                    {% endfor %}
                  </div>
                </div>
                <button type="button" class="btn btn-primary" id="apply-btn">Apply</button>
              </form>
            </div>
          </div>
        </div>
        <div class="col-md-9" id="main-panel">
          <div class="card">
            <div class="card-header">Results</div>
            <div class="card-body">
              <div id="results-area">
                <p class="text-muted">Use the controls to select a sample/contig and modules, then press Apply.</p>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>

    <script src="https://code.jquery.com/jquery-3.5.1.slim.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/axios/dist/axios.min.js"></script>
    <script>
      document.getElementById('toggle-btn').addEventListener('click', function() {
        var controls = document.getElementById('controls');
        if (controls.style.display === 'none') {
          controls.style.display = 'block';
          this.textContent = 'Hide';
        } else {
          controls.style.display = 'none';
          this.textContent = 'Show';
        }
      });

      document.getElementById('apply-btn').addEventListener('click', function() {
        var sample = document.getElementById('sample-select').value;
        var contig = document.getElementById('contig-select').value;
        var modules = [];
        document.querySelectorAll('.module-checkbox:checked').forEach(function(cb){ modules.push(cb.value); });

        var payload = { sample: sample, contig: contig, modules: modules };
        document.getElementById('results-area').innerHTML = '<p>Loading...</p>';

        axios.post('/fetch_features', payload).then(function(resp){
          var data = resp.data;
          // For now show a pretty-printed JSON extract in the main panel
          var pretty = JSON.stringify(data, null, 2);
          document.getElementById('results-area').innerHTML = '<pre class="json">' + escapeHtml(pretty) + '</pre>';
        }).catch(function(err){
          console.error(err);
          document.getElementById('results-area').innerHTML = '<div class="alert alert-danger">Error fetching features (see console)</div>';
        });
      });

      function escapeHtml(unsafe) {
        return unsafe
          .replace(/&/g, "&amp;")
          .replace(/</g, "&lt;")
          .replace(/>/g, "&gt;")
          .replace(/\"/g, "&quot;")
          .replace(/'/g, "&#039;");
      }
    </script>
  </body>
</html>
"""


def query_db(query, params=()):
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()
    cur.execute(query, params)
    rows = cur.fetchall()
    conn.close()
    return rows


@app.route('/')
def index():
    # populate sample and contig lists and module list from DB
    samples = [r['Sample_name'] for r in query_db('SELECT Sample_name FROM Sample ORDER BY Sample_name')]
    contigs = [r['Contig_name'] for r in query_db('SELECT Contig_name FROM Contig ORDER BY Contig_name')]
    modules = [r[0] for r in query_db('SELECT DISTINCT Module FROM Variable WHERE Module IS NOT NULL ORDER BY Module')]
    return render_template_string(HTML_TEMPLATE, samples=samples, contigs=contigs, modules=modules)


@app.route('/fetch_features', methods=['POST'])
def fetch_features():
    data = request.get_json()
    sample = data.get('sample')
    contig = data.get('contig')
    modules = data.get('modules') or []

    if not sample or not contig or not modules:
        return jsonify({'error': 'sample, contig and at least one module required'}), 400

    # open connection once per request
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    cur = conn.cursor()

    # map sample / contig to ids
    cur.execute('SELECT Sample_id FROM Sample WHERE Sample_name=?', (sample,))
    srow = cur.fetchone()
    if not srow:
        conn.close()
        return jsonify({'error': f'Sample {sample} not found'}), 404
    sample_id = srow['Sample_id']

    cur.execute('SELECT Contig_id FROM Contig WHERE Contig_name=?', (contig,))
    crow = cur.fetchone()
    if not crow:
        conn.close()
        return jsonify({'error': f'Contig {contig} not found'}), 404
    contig_id = crow['Contig_id']

    result = {}

    # For each module, find variables and fetch a small extract of data for the sample/contig
    for module in modules:
        cur.execute('SELECT Variable_name, Feature_table_name FROM Variable WHERE Module=?', (module,))
        vars_rows = cur.fetchall()
        module_data = {}
        for vr in vars_rows:
            var_name = vr['Variable_name']
            feature_table = vr['Feature_table_name']
            # Basic safety: allow only alnum + underscore in table names
            if not feature_table.replace('_', '').isalnum():
                continue
            # Query a small extract (limit 200 rows)
            try:
                q = f"SELECT Position, Value FROM {feature_table} WHERE Sample_id=? AND Contig_id=? ORDER BY Position LIMIT 200"
                cur.execute(q, (sample_id, contig_id))
                rows = [{'position': int(r['Position']), 'value': float(r['Value'])} for r in cur.fetchall()]
            except Exception as e:
                rows = [{'error': str(e)}]
            module_data[var_name] = rows
        result[module] = module_data

    conn.close()
    # For now return the raw data extract to the client
    return jsonify(result)


def main():
    global DB_PATH
    parser = argparse.ArgumentParser(description='Serve a small HTML front-end to inspect feature database')
    parser.add_argument('--db', required=True, help='Path to sqlite DB produced by generate_database.py')
    parser.add_argument('--outfile', required=False, help='If set, render the HTML to this file and exit (no server)')
    parser.add_argument('--host', default='127.0.0.1', help='Host to bind to')
    parser.add_argument('--port', type=int, default=5000, help='Port to bind to')
    args = parser.parse_args()

    DB_PATH = args.db

    # sanity check DB exists
    if not os.path.exists(DB_PATH):
        print(f"Database not found: {DB_PATH}. Run generate_database.py <db> first.")
        raise SystemExit(1)

    # If an outfile was requested, render template with DB values and write static HTML then exit.
    if args.outfile:
        samples = [r['Sample_name'] for r in query_db('SELECT Sample_name FROM Sample ORDER BY Sample_name')]
        contigs = [r['Contig_name'] for r in query_db('SELECT Contig_name FROM Contig ORDER BY Contig_name')]
        modules = [r[0] for r in query_db('SELECT DISTINCT Module FROM Variable WHERE Module IS NOT NULL ORDER BY Module')]
        rendered = render_template_string(HTML_TEMPLATE, samples=samples, contigs=contigs, modules=modules)
        with open(args.outfile, 'w', encoding='utf-8') as fh:
            fh.write(rendered)
        print(f'Wrote static HTML to: {args.outfile}')
        return

    # Start the Flask app
    app.run(host=args.host, port=args.port, debug=True)


if __name__ == '__main__':
    main()
