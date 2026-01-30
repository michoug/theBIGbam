import param
from panel.custom import JSComponent


class SearchableSelect(JSComponent):
    """Searchable select dropdown using Tom Select."""

    value = param.String(default="")
    options = param.List(default=[])
    placeholder = param.String(default="")

    _stylesheets = [
        "https://cdn.jsdelivr.net/npm/tom-select@2.4.1/dist/css/tom-select.css",
        # Override styles to match Bokeh widget look
        """
        .ts-wrapper { width: 100%; }
        .ts-wrapper .ts-control { border: 1px solid #ccc; border-radius: 4px; padding: 4px 8px; min-height: 31px; font-size: 14px; }
        """
    ]

    _esm = """
    import TomSelect from "https://esm.sh/tom-select@2.4.1";

    export function render({ model }) {
        const container = document.createElement('div');
        container.style.width = '100%';
        container.style.margin = '5px 10px 10px 0';
        const select = document.createElement('select');
        select.setAttribute('placeholder', model.placeholder);
        container.appendChild(select);

        const buildOptions = () => {
            return model.options.map(o => ({value: o, text: o}));
        };

        const ts = new TomSelect(select, {
            create: false,
            maxOptions: null,
            placeholder: model.placeholder,
            options: buildOptions(),
            items: model.value ? [model.value] : [],
            onChange: (val) => { model.value = val; }
        });

        model.on('options', () => {
            const currentVal = model.value;
            ts.clearOptions();
            ts.addOptions(buildOptions());
            // Restore value if still valid, otherwise clear
            if (model.options.includes(currentVal)) {
                ts.setValue(currentVal, true);
            } else {
                ts.setValue('', true);
            }
        });

        model.on('value', () => {
            ts.setValue(model.value, true);
        });

        return container;
    }
    """
