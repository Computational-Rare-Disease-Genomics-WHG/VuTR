"""WSGI to run the application"""

from app import create_app  # pylint: disable=C0114

app = create_app()

if __name__ == '__main__':
    app.run()
