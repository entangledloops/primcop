from dash import Dash

from src.components.layout import create_layout


def main() -> None:
    app = Dash()
    app.title = "Potential Prion Sequence Analyzer"
    app.layout = create_layout(app)
    app.run_server()


if __name__ == "__main__":
    main()