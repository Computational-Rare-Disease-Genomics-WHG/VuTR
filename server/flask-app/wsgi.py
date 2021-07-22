from app import * 

app = create_app("prod")

if __name__ == "__main__": 
    app.run()
